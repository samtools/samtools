#include <time.h>
#include <stdio.h>
#include <netdb.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include "knetfile.h"

static int socket_wait(int fd, int is_read)
{
	fd_set fds, *fdr = 0, *fdw = 0;
	struct timeval tv;
	int ret;
	tv.tv_sec = 5; tv.tv_usec = 0; // 5 seconds time out
	FD_ZERO(&fds);
	FD_SET(fd, &fds);
	if (is_read) fdr = &fds;
	else fdw = &fds;
	ret = select(fd+1, fdr, fdw, 0, &tv);
	if (ret == -1) perror("select");
	return ret;
}

static int socket_connect(const char *host, const char *port)
{
#define __err_connect(func) do { perror(func); freeaddrinfo(res); return -1; } while (0)

	int on = 1, fd;
	struct linger lng = { 0, 0 };
	struct addrinfo hints, *res;
	memset(&hints, 0, sizeof(struct addrinfo));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	if (getaddrinfo(host, port, &hints, &res) != 0) __err_connect("getaddrinfo");
	if ((fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)) == -1) __err_connect("socket");
	if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on)) == -1) __err_connect("setsockopt");
	if (setsockopt(fd, SOL_SOCKET, SO_LINGER, &lng, sizeof(lng)) == -1) __err_connect("setsockopt");
	if (connect(fd, res->ai_addr, res->ai_addrlen) != 0) __err_connect("connect");
	freeaddrinfo(res);
	return fd;
}

static off_t my_read(int fd, void *buf, off_t len)
{
	off_t rest = len, curr, l = 0;
	while (rest) {
		if (socket_wait(fd, 1) <= 0) break; // socket is not ready for reading
		curr = read(fd, buf + l, rest);
		if (curr == 0) break;
		l += curr; rest -= curr;
	}
	return l;
}

/*************************
 * FTP specific routines *
 *************************/

static int kftp_get_response(knetFile *ftp)
{
	unsigned char c;
	int n = 0;
	char *p;
	if (socket_wait(ftp->ctrl_fd, 1) <= 0) return 0;
	while (read(ftp->ctrl_fd, &c, 1)) { // FIXME: this is *VERY BAD* for unbuffered I/O
		//fputc(c, stderr);
		if (n >= ftp->max_response) {
			ftp->max_response = ftp->max_response? ftp->max_response<<1 : 256;
			ftp->response = realloc(ftp->response, ftp->max_response);
		}
		ftp->response[n++] = c;
		if (c == '\n') {
			if (n >= 4 && isdigit(ftp->response[0]) && isdigit(ftp->response[1]) && isdigit(ftp->response[2])
				&& ftp->response[3] != '-') break;
			n = 0;
			continue;
		}
	}
	if (n < 2) return -1;
	ftp->response[n-2] = 0;
	return strtol(ftp->response, &p, 0);
}

static int kftp_send_cmd(knetFile *ftp, const char *cmd, int is_get)
{
	if (socket_wait(ftp->ctrl_fd, 0) <= 0) return -1; // socket is not ready for writing
	write(ftp->ctrl_fd, cmd, strlen(cmd));
	return is_get? kftp_get_response(ftp) : 0;
}

static int kftp_pasv_prep(knetFile *ftp)
{
	char *p;
	int v[6];
	kftp_send_cmd(ftp, "PASV\r\n", 1);
	for (p = ftp->response; *p && *p != '('; ++p);
	if (*p != '(') return -1;
	++p;
	sscanf(p, "%d,%d,%d,%d,%d,%d", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]);
	memcpy(ftp->pasv_ip, v, 4 * sizeof(int));
	ftp->pasv_port = (v[4]<<8&0xff00) + v[5];
	return 0;
}


static int kftp_pasv_connect(knetFile *ftp)
{
	char host[80], port[10];
	if (ftp->pasv_port == 0) {
		fprintf(stderr, "[kftp_pasv_connect] kftp_pasv_prep() is not called before hand.\n");
		return -1;
	}
	sprintf(host, "%d.%d.%d.%d", ftp->pasv_ip[0], ftp->pasv_ip[1], ftp->pasv_ip[2], ftp->pasv_ip[3]);
	sprintf(port, "%d", ftp->pasv_port);
	ftp->fd = socket_connect(host, port);
	if (ftp->fd == -1) return -1;
	return 0;
}

int kftp_connect(knetFile *ftp)
{
	ftp->ctrl_fd = socket_connect(ftp->host, ftp->port);
	if (ftp->ctrl_fd == -1) return -1;
	kftp_get_response(ftp);
	kftp_send_cmd(ftp, "USER anonymous\r\n", 1);
	kftp_send_cmd(ftp, "PASS kftp@\r\n", 1);
	kftp_send_cmd(ftp, "TYPE I\r\n", 1);
	return 0;
}

int kftp_reconnect(knetFile *ftp)
{
	if (ftp->ctrl_fd >= 0) {
		close(ftp->ctrl_fd);
		ftp->ctrl_fd = -1;
	}
	close(ftp->fd);
	return kftp_connect(ftp);
}

// initialize ->type, ->host and ->retr
knetFile *kftp_parse_url(const char *fn, const char *mode)
{
	knetFile *fp;
	char *p;
	int l;
	if (strstr(fn, "ftp://") != fn) return 0;
	for (p = (char*)fn + 6; *p && *p != '/'; ++p);
	if (*p != '/') return 0;
	l = p - fn - 6;
	fp = calloc(1, sizeof(knetFile));
	fp->type = KNF_TYPE_FTP;
	fp->fd = -1;
	fp->port = strdup("ftp");
	fp->host = calloc(l + 1, 1);
	if (strchr(mode, 'c')) fp->no_reconnect = 1;
	strncpy(fp->host, fn + 6, l);
	fp->retr = calloc(strlen(p) + 8, 1);
	sprintf(fp->retr, "RETR %s\r\n", p);
	fp->seek_offset = -1;
	return fp;
}
// place ->fd at offset off
int kftp_connect_file(knetFile *fp)
{
	int ret;
	if (fp->fd >= 0) {
		close(fp->fd);
		if (fp->no_reconnect) kftp_get_response(fp);
	}
	kftp_pasv_prep(fp);
	if (fp->offset) {
		char tmp[32];
		sprintf(tmp, "REST %lld\r\n", (long long)fp->offset);
		kftp_send_cmd(fp, tmp, 1);
	}
	kftp_send_cmd(fp, fp->retr, 0);
	kftp_pasv_connect(fp);
	ret = kftp_get_response(fp);
	if (ret != 150) {
		fprintf(stderr, "[kftp_connect_file] %s\n", fp->response);
		close(fp->fd);
		fp->fd = -1;
		return -1;
	}
	fp->is_ready = 1;
	return 0;
}

/**************************
 * HTTP specific routines *
 **************************/

knetFile *khttp_parse_url(const char *fn, const char *mode)
{
	knetFile *fp;
	char *p, *proxy, *q;
	int l;
	if (strstr(fn, "http://") != fn) return 0;
	// set ->http_host
	for (p = (char*)fn + 7; *p && *p != '/'; ++p);
	l = p - fn - 7;
	fp = calloc(1, sizeof(knetFile));
	fp->http_host = calloc(l + 1, 1);
	strncpy(fp->http_host, fn + 7, l);
	fp->http_host[l] = 0;
	for (q = fp->http_host; *q && *q != ':'; ++q);
	if (*q == ':') *q++ = 0;
	// get http_proxy
	proxy = getenv("http_proxy");
	// set ->host, ->port and ->path
	if (proxy == 0) {
		fp->host = strdup(fp->http_host); // when there is no proxy, server name is identical to http_host name.
		fp->port = strdup(*q? q : "http");
		fp->path = strdup(*p? p : "/");
	} else {
		fp->host = (strstr(proxy, "http://") == proxy)? strdup(proxy + 7) : strdup(proxy);
		for (q = fp->host; *q && *q != ':'; ++q);
		if (*q == ':') *q++ = 0; 
		fp->port = strdup(*q? q : "http");
		fp->path = strdup(fn);
	}
	fp->type = KNF_TYPE_HTTP;
	fp->ctrl_fd = fp->fd = -1;
	fp->seek_offset = -1;
	return fp;
}

int khttp_connect_file(knetFile *fp)
{
	int ret, l = 0;
	char *buf, *p;
	if (fp->fd >= 0) close(fp->fd);
	fp->fd = socket_connect(fp->host, fp->port);
	buf = calloc(0x10000, 1); // FIXME: I am lazy... But in principle, 64KB should be large enough.
	l += sprintf(buf + l, "GET %s HTTP/1.0\r\nHost: %s\r\n", fp->path, fp->http_host);
	if (fp->offset)
		l += sprintf(buf + l, "Range: bytes=%lld-\r\n", (long long)fp->offset);
	l += sprintf(buf + l, "\r\n");
	write(fp->fd, buf, l);
	l = 0;
	while (read(fp->fd, buf + l, 1)) { // read HTTP header; FIXME: bad efficiency
		if (buf[l] == '\n' && l >= 3)
			if (strncmp(buf + l - 3, "\r\n\r\n", 4) == 0) break;
		++l;
	}
	buf[l] = 0;
	if (l < 14) { // prematured header
		close(fp->fd);
		fp->fd = -1;
		return -1;
	}
	ret = strtol(buf + 8, &p, 0); // HTTP return code
	if (ret == 200 && fp->offset) { // 200 (complete result); then skip beginning of the file
		off_t rest = fp->offset;
		while (rest) {
			off_t l = rest < 0x10000? rest : 0x10000;
			rest -= my_read(fp->fd, buf, l);
		}
	} else if (ret != 206 && ret != 200) {
		free(buf);
		fprintf(stderr, "[khttp_connect_file] fail to open file (HTTP code: %d).\n", ret);
		close(fp->fd);
		fp->fd = -1;
		return -1;
	}
	free(buf);
	fp->is_ready = 1;
	return 0;
}

/********************
 * Generic routines *
 ********************/

knetFile *knet_open(const char *fn, const char *mode)
{
	knetFile *fp = 0;
	if (mode[0] != 'r') {
		fprintf(stderr, "[kftp_open] only mode \"r\" is supported.\n");
		return 0;
	}
	if (strstr(fn, "ftp://") == fn) {
		fp = kftp_parse_url(fn, mode);
		if (fp == 0) return 0;
		if (kftp_connect(fp) == -1) {
			knet_close(fp);
			return 0;
		}
		kftp_connect_file(fp);
	} else if (strstr(fn, "http://") == fn) {
		fp = khttp_parse_url(fn, mode);
		if (fp == 0) return 0;
		khttp_connect_file(fp);
	} else { // local file
		int fd = open(fn, O_RDONLY);
		if (fd == -1) {
			perror("open");
			return 0;
		}
		fp = (knetFile*)calloc(1, sizeof(knetFile));
		fp->type = KNF_TYPE_LOCAL;
		fp->fd = fd;
		fp->ctrl_fd = -1;
	}
	if (fp && fp->fd < 0) {
		knet_close(fp);
		return 0;
	}
	return fp;
}

knetFile *knet_dopen(int fd, const char *mode)
{
	knetFile *fp = (knetFile*)calloc(1, sizeof(knetFile));
	fp->type = KNF_TYPE_LOCAL;
	fp->fd = fd;
	return fp;
}

off_t knet_read(knetFile *fp, void *buf, off_t len)
{
	off_t l = 0;
	if (fp->fd < 0) return 0;
	if (fp->type == KNF_TYPE_FTP) {
		if (fp->is_ready == 0) {
			if (!fp->no_reconnect) kftp_reconnect(fp);
			kftp_connect_file(fp);
		}
	} else if (fp->type == KNF_TYPE_HTTP) {
		if (fp->is_ready == 0)
			khttp_connect_file(fp);
	}
	l = my_read(fp->fd, buf, len);
	fp->offset += l;
	return l;
}

int knet_seek(knetFile *fp, off_t off, int whence)
{
	if (whence == SEEK_SET && off == fp->offset) return 0;
	if (fp->type == KNF_TYPE_LOCAL) {
		off_t offset = lseek(fp->fd, off, whence);
		if (offset == -1) {
			perror("lseek");
			return -1;
		}
		fp->offset = offset;
		return 0;
	} else if (fp->type == KNF_TYPE_FTP || fp->type == KNF_TYPE_HTTP) {
		if (whence != SEEK_SET) { // FIXME: we can surely allow SEEK_CUR and SEEK_END in future
			fprintf(stderr, "[knet_seek] only SEEK_SET is supported for FTP/HTTP. Offset is unchanged.\n");
			return -1;
		}
		fp->offset = off;
		fp->is_ready = 0;
		return 0;
	}
	return -1;
}

int knet_close(knetFile *fp)
{
	if (fp == 0) return 0;
	if (fp->ctrl_fd >= 0) close(fp->ctrl_fd); // FTP specific
	if (fp->fd >= 0) close(fp->fd);
	free(fp->host); free(fp->port);
	free(fp->response); free(fp->retr); // FTP specific
	free(fp->path); free(fp->http_host); // HTTP specific
	free(fp);
	return 0;
}

#ifdef KNETFILE_MAIN
int main(void)
{
	char *buf;
	knetFile *fp;
	int type = 4, l;
	buf = calloc(0x100000, 1);
	if (type == 0) {
		fp = knet_open("knetfile.c", "r");
		knet_seek(fp, 1000, SEEK_SET);
	} else if (type == 1) { // NCBI FTP, large file
		fp = knet_open("ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12878/alignment/NA12878.chrom6.SLX.SRP000032.2009_06.bam", "r");
		knet_seek(fp, 2500000000ll, SEEK_SET);
		l = knet_read(fp, buf, 255);
	} else if (type == 2) {
		fp = knet_open("ftp://ftp.sanger.ac.uk/pub4/treefam/tmp/index.shtml", "r");
		knet_seek(fp, 1000, SEEK_SET);
	} else if (type == 3) {
		fp = knet_open("http://www.sanger.ac.uk/Users/lh3/index.shtml", "r");
		knet_seek(fp, 1000, SEEK_SET);
	} else if (type == 4) {
		fp = knet_open("http://www.sanger.ac.uk/Users/lh3/ex1.bam", "r");
		knet_read(fp, buf, 10000);
		knet_seek(fp, 20000, SEEK_SET);
		knet_seek(fp, 10000, SEEK_SET);
		l = knet_read(fp, buf+10000, 10000000) + 10000;
	}
	if (type != 4 && type != 1) {
		knet_read(fp, buf, 255);
		buf[255] = 0;
		printf("%s\n", buf);
	} else write(fileno(stdout), buf, l);
	knet_close(fp);
	free(buf);
	return 0;
}
#endif
