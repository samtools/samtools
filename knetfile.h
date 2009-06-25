#ifndef KNETFILE_H
#define KNETFILE_H

#include <stdint.h>
#include <fcntl.h>

// FIXME: currently I/O is unbuffered

#define KNF_TYPE_LOCAL 1
#define KNF_TYPE_FTP   2
#define KNF_TYPE_HTTP  3

typedef struct knetFile_s {
	int type, fd;
	int64_t offset;
	char *host;

	// the following are for FTP only
	int ctrl_fd, pasv_ip[4], pasv_port, max_response, no_reconnect;
	char *response, *retr;
} knetFile;

#define knet_tell(fp) ((fp)->offset)
#define knet_fileno(fp) ((fp)->fd)

#ifdef __cplusplus
extern "C" {
#endif

	knetFile *knet_open(const char *fn, const char *mode);
	knetFile *knet_dopen(int fd, const char *mode);
	off_t knet_read(knetFile *fp, void *buf, off_t len);
	off_t knet_seek(knetFile *fp, off_t off, int whence);
	int knet_close(knetFile *fp);

#ifdef __cplusplus
}
#endif

#endif
