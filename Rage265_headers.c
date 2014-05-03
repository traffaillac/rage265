/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */

/**
 * Order of inclusion matters, lower files use functions and variables from
 * higher files. Compiling all dependencies along with the main library has
 * advantages not met by Whole Program Optimisation:
 * _ simplicity: no need to expose internal APIs between compilation units,
 *   avoiding the .c/.h couples thus reducing the number of files;
 * _ transparency: compilation needs not be hidden behind a Makefile, and the
 *   simple command used to build the library can be tuned as will;
 * _ output control: static functions do not appear in the resulting archive,
 *   without having to strip them afterwards.
 */
#include "Rage265_common.h"



/**
 * Find the next 00 00 0n pattern, pointing to its first byte or the buffer end.
 */
#ifdef __SSSE3__
const uint8_t *Rage265_find_start_code(const uint8_t *cur, size_t len, unsigned int n) {
	const uint8_t *end = cur + len;
	while (cur < end) {
		/* Perform a per-byte search until cur is aligned. */
		const __m128i *x = (__m128i *)((uintptr_t)cur & -sizeof(*x)) + 1;
		const uint8_t *stop = ((uint8_t *)x + 2 < end) ? (uint8_t *)x + 2 : end;
		for (unsigned int start_code = -1; cur < stop; cur++) {
			start_code = ((start_code & 0xffff) << 8) | *cur;
			if (start_code == n)
				return cur - 2;
		}
		
		/* Skip words without a zero odd byte. */
		while (x < (__m128i *)((uintptr_t)end & -sizeof(*x)) && (_mm_movemask_epi8(_mm_cmpeq_epi8(*x, _mm_setzero_si128())) & 0xaaaa) == 0)
			cur = (uint8_t *)++x;
	}
	return cur;
}
#endif



/**
 * Parse a NAL unit, returning a bitfield of error codes.
 */
unsigned int Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *cur, size_t len) {
	unsigned int nal_unit_header = (cur[0] << 8) | cur[1];
	
}
