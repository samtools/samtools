#include <errno.h>

void xfreopen(const char *path, const char *mode, FILE *stream)
{
	if (freopen(path, mode, stream) == NULL) {
		fprintf(stderr, __FILE__": error reopening %s: %s\n",
				path, strerror(errno));
		exit(2);
	}
}

