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

static int kftp_get_response(knetFile *ftp)
{
	unsigned char c;
	int n = 0;
	char *p;
	while (read(ftp->ctrl_fd, &c, 1)) { // FIXME: this is *VERY BAD* for unbuffered I/O
//		fputc(c, stderr);
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
#define __err_pasv_connect(func) do { perror(func); freeaddrinfo(res); return -1; } while (0)

	struct addrinfo hints, *res;
	struct linger lng = { 0, 0 };
	int on = 1;
	char host[80], port[10];

	if (ftp->pasv_port == 0) {
		fprintf(stderr, "[kftp_pasv_connect] kftp_pasv_prep() is not called before hand.\n");
		return -1;
	}
	memset(&hints, 0, sizeof(struct addrinfo));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	sprintf(host, "%d.%d.%d.%d", ftp->pasv_ip[0], ftp->pasv_ip[1], ftp->pasv_ip[2], ftp->pasv_ip[3]);
	sprintf(port, "%d", ftp->pasv_port);
	if (getaddrinfo(host, port, &hints, &res) != 0) { perror("getaddrinfo"); return -1; }
	if ((ftp->fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)) == -1) __err_pasv_connect("socket");
	if (setsockopt(ftp->fd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on)) == -1) __err_pasv_connect("setsockopt");
	if (setsockopt(ftp->fd, SOL_SOCKET, SO_LINGER, &lng, sizeof(lng)) == -1) __err_pasv_connect("setsockopt");
	if (connect(ftp->fd, res->ai_addr, res->ai_addrlen) != 0) __err_pasv_connect("connect");
	freeaddrinfo(res);
	return 0;
}

int kftp_connect(knetFile *ftp)
{
#define __err_connect(func) do { perror(func); return -1; } while (0)

	int on = 1;
	{ // open socket
		struct addrinfo hints, *res;
		memset(&hints, 0, sizeof(struct addrinfo));
		hints.ai_family = AF_UNSPEC;
		hints.ai_socktype = SOCK_STREAM;
		if (getaddrinfo(ftp->host, "21", &hints, &res) != 0) __err_connect("getaddrinfo");
		if ((ftp->ctrl_fd = socket(res->ai_family, res->ai_socktype, res->ai_protocol)) == -1) __err_connect("socket");
		if (setsockopt(ftp->ctrl_fd, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on)) == -1) __err_connect("setsockopt");
		if (connect(ftp->ctrl_fd, res->ai_addr, res->ai_addrlen) != 0) __err_connect("connect");
		freeaddrinfo(res);
		kftp_get_response(ftp);
	}
	{ // login
		kftp_send_cmd(ftp, "USER anonymous\r\n", 1);
		kftp_send_cmd(ftp, "PASS kftp@\r\n", 1);
		kftp_send_cmd(ftp, "TYPE I\r\n", 1);
	}
	return 0;
}

int kftp_reconnect(knetFile *ftp)
{
	if (ftp->ctrl_fd) {
		close(ftp->ctrl_fd);
		ftp->ctrl_fd = 0;
	}
	close(ftp->fd);
	return kftp_connect(ftp);
}

// initialize ->type, ->host and ->retr
knetFile *kftp_prep(const char *fn, const char *mode)
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
	if (fp->fd) {
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
	kftp_get_response(fp);
	fp->is_ready = 1;
	return 0;
}

knetFile *knet_open(const char *fn, const char *mode)
{
	knetFile *fp = 0;
	if (mode[0] != 'r') {
		fprintf(stderr, "[kftp_open] only mode \"r\" is supported.\n");
		return 0;
	}
	if (strstr(fn, "ftp://") == fn) {
		fp = kftp_prep(fn, mode);
		if (fp == 0) return 0;
		if (kftp_connect(fp) == -1) {
			knet_close(fp);
			return 0;
		}
		kftp_connect_file(fp);
	} else {
		int fd = open(fn, O_RDONLY);
		if (fd == -1) {
			perror("open");
			return 0;
		}
		fp = (knetFile*)calloc(1, sizeof(knetFile));
		fp->type = KNF_TYPE_LOCAL;
		fp->fd = fd;
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
	if (fp->type == KNF_TYPE_LOCAL) {
		l = read(fp->fd, buf, len);
		fp->offset += l;
	} else {
		off_t rest = len, curr;
		if (fp->is_ready == 0) {
			if (!fp->no_reconnect) kftp_reconnect(fp);
			kftp_connect_file(fp);
			fp->is_ready = 1;
		}
		while (rest) {
			curr = read(fp->fd, buf + l, rest);
			if (curr == 0) break; // FIXME: end of file or bad network? I do not know...
			l += curr; rest -= curr;
		}
		fp->offset += l;
	}
	return l;
}

int knet_seek(knetFile *fp, off_t off, int whence)
{
	if (fp->type == KNF_TYPE_LOCAL) {
		if (lseek(fp->fd, off, whence) == -1) {
			perror("lseek");
			return -1;
		}
		fp->offset = off;
		return 0;
	}
	if (fp->type == KNF_TYPE_FTP) {
		if (whence != SEEK_SET) { // FIXME: we can surely allow SEEK_CUR and SEEK_END in future
			fprintf(stderr, "[knet_seek] only SEEK_SET is supported for FTP. Offset is unchanged.\n");
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
	if (fp->ctrl_fd > 0) close(fp->ctrl_fd);
	if (fp->fd > 0) close(fp->fd);
	free(fp->response); free(fp->retr); free(fp->host);
	free(fp);
	return 0;
}

#ifdef KNETFILE_MAIN
int main(void)
{
	char buf[256];
	knetFile *fp;
//	fp = knet_open("ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12878/alignment/NA12878.chrom6.SLX.SRP000032.2009_06.bam", "r"); knet_seek(fp, 2500000000ll, SEEK_SET);
	fp = knet_open("ftp://ftp.sanger.ac.uk/pub4/treefam/tmp/index.shtml", "r"); knet_seek(fp, 2000, SEEK_SET);
//	fp = knet_open("knetfile.c", "r"); knet_seek(fp, 2000, SEEK_SET);
	knet_read(fp, buf, 255);
	buf[255] = 0;
	printf("%s\n", buf);
	knet_close(fp);
	return 0;
}
#endif
