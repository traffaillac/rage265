#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#define DEBUG 2
#include "Rage265.c"

int main() {
	/* Memory-map the whole file. */
	struct stat st;
	fstat(0, &st);
	const uint8_t *file = mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, 0, 0);
	assert(file!=MAP_FAILED);
	setvbuf(stdout, NULL, _IONBF, BUFSIZ);
	
	/* Parse and dump the file to HTML. */
	printf("<!doctype html>\n"
		"<html>\n"
		"<head><meta charset=\"UTF-8\"/><title>Rage265 dump</title></head>\n"
		"<body>\n");
	Rage265_ctx r = {0};
	for (size_t len, i = 4; i < st.st_size; i += len + 3) {
		len = Rage265_find_start_code(file + i, st.st_size - i, 1);
		const Rage265_picture *p = Rage265_parse_NAL(&r, file + i, len);
		if (p != NULL)
			printf("<p>Output picture %d</p>\n", p->PicOrderCntVal);
	}
	printf("</body>\n"
		"</html>\n");
	return 0;
}
