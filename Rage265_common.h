/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */
#ifndef RAGE265_COMMON_H
#define RAGE265_COMMON_H

#if DEBUG >= 1
#include <stdio.h>
static inline const char *red_if(int cond) { return (cond) ? " style=\"color: red\"" : ""; }
#else
#define printf(...) ((void)0)
#define NDEBUG 1
#endif
#if DEBUG != 2
#define fprintf(...) ((void)0)
#endif

#include <assert.h>
#include <endian.h>
#include <limits.h>

#include "Rage265.h"

#ifndef WORD_BIT
#if INT_MAX == 2147483647
#define WORD_BIT 32
#endif
#endif
#ifndef LONG_BIT
#if LONG_MAX == 2147483647
#define LONG_BIT 32
#elif LONG_MAX == 9223372036854775807
#define LONG_BIT 64
#endif
#endif

#ifdef __SSSE3__
#define _mm_movpi64_pi64 _mm_movpi64_epi64
#ifdef __SSE4_1__
#include <smmintrin.h>
#else
#include <tmmintrin.h>
#define _mm_extract_epi32(a, i) \
	_mm_cvtsi128_si32(_mm_shuffle_epi32(a, _MM_SHUFFLE(i, i, i, i)))
static inline __m128i _mm_packus_epi32(__m128i a, __m128i b) {
	return _mm_max_epi16(_mm_packs_epi32(a, b), _mm_setzero_si128());
}
static inline __m128i _mm_mullo_epi32(__m128i a, __m128i b) {
	__m128i x0 = _mm_shuffle_epi32(a, _MM_SHUFFLE(0, 3, 0, 1));
	__m128i x1 = _mm_shuffle_epi32(b, _MM_SHUFFLE(0, 3, 0, 1));
	__m128i x2 = _mm_mul_epu32(a, b);
	__m128i x3 = _mm_mul_epu32(x0, x1);
	__m128 x4 = _mm_shuffle_ps((__m128)x2, (__m128)x3, _MM_SHUFFLE(2, 0, 2, 0));
	return _mm_shuffle_epi32((__m128i)x4, _MM_SHUFFLE(3, 1, 2, 0));
}
#endif
#endif

static inline long min(long a, long b) { return (a < b) ? a : b; }
static inline long max(long a, long b) { return (a > b) ? a : b; }



/**
 * 9.2 - Exp-Golomb parsing
 */
#if ULLONG_MAX == 18446744073709551615U
#define clz64 __builtin_clzll
#endif

static inline __attribute__((always_inline)) unsigned int get_ue8(const uint8_t * restrict CPB, unsigned int * restrict shift) {
	uint16_t buf = (((CPB[*shift / 8] << 8) | CPB[*shift / 8 + 1]) << (*shift % 8));
	unsigned int leadingZeroBits = __builtin_clz(buf | 0x0400) - WORD_BIT + 16;
	*shift += 2 * leadingZeroBits + 1;
	return (buf >> (16 - (2 * leadingZeroBits + 1))) - 1;
}

static inline __attribute__((always_inline)) unsigned int get_ue32(const uint8_t * restrict CPB, unsigned int * restrict shift) {
	unsigned int msb = htobe32(((uint32_t *)CPB)[*shift / 32]);
	unsigned int lsb = htobe32(((uint32_t *)CPB)[(*shift + 31) / 32]);
	uint32_t buf = (msb << (*shift % 32)) | (lsb >> (-*shift % 32));
	unsigned int leadingZeroBits = __builtin_clz(buf | 0x00010000) - WORD_BIT + 32;
	*shift += 2 * leadingZeroBits + 1;
	return (buf >> (32 - (2 * leadingZeroBits + 1))) - 1;
}

static unsigned int get_ue64(const uint8_t * restrict CPB, unsigned int * restrict shift) {
	uint64_t msb = htobe64(((uint64_t *)CPB)[*shift / 64]);
	uint64_t lsb = htobe64(((uint64_t *)CPB)[(*shift + 63) / 64]);
	uint64_t buf = (msb << (*shift % 64)) | (lsb >> (-*shift % 64));
	unsigned int leadingZeroBits = clz64(buf | 0x0000000100000000);
	*shift += 2 * leadingZeroBits + 1;
	return (buf >> (64 - (2 * leadingZeroBits + 1))) - 1;
}

static inline __attribute__((always_inline)) unsigned int get_ue(const uint8_t * restrict CPB, unsigned int * restrict shift, unsigned int upper) {
	assert(upper<4294967295);
	unsigned int res = (upper <= 31) ? get_ue8(CPB, shift) : (upper <= 65534) ? get_ue32(CPB, shift) : get_ue64(CPB, shift);
	return (res < upper) ? res : upper;
}

static inline __attribute__((always_inline)) int get_se(const uint8_t * restrict CPB, unsigned int * restrict shift, int lower, int upper) {
	unsigned int codeNum = get_ue(CPB, shift, max(-lower * 2, upper * 2 - 1));
	unsigned int abs = (codeNum+ 1) / 2;
	unsigned int sign = (codeNum % 2) - 1;
	int res = (abs ^ sign) - sign;
	return (-lower * 2 < upper * 2 - 1) ? max(res, lower) : (-lower * 2 > upper * 2 - 1) ? min(res, upper) : res;
}

static inline __attribute__((always_inline)) unsigned int get_uv(const uint8_t * restrict CPB, unsigned int * restrict shift, unsigned int v) {
	unsigned int msb = htobe32(((uint32_t *)CPB)[*shift / 32]);
	unsigned int lsb = htobe32(((uint32_t *)CPB)[(*shift + 31) / 32]);
	uint32_t buf = (msb << (*shift % 32)) | (lsb >> (-*shift % 32));
	*shift += v;
	return buf >> (32 - v);
}

static inline __attribute__((always_inline)) unsigned int get_u1(const uint8_t * restrict CPB, unsigned int * restrict shift) {
	uint8_t buf = CPB[*shift / 8];
	return (buf >> (7 - *shift++ % 8)) & 1;
}



/**
 * 9.3 - CABAC parsing
 */
typedef struct {
	unsigned long ivlCurrRange;
	unsigned long ivlOffset;
	const uint8_t *CPB;
	unsigned int shift;
	unsigned int lim;
} CABAC_ctx;

#if ULONG_MAX == 4294967295U
#define htobe htobe32
#elif ULONG_MAX == 18446744073709551615U
#define htobe htobe64
#endif

static inline void renorm(CABAC_ctx *c, unsigned int v) {
	unsigned long buf = 0;
	if (c->shift < c->lim) {
		unsigned long msb = htobe(((unsigned long *)c->CPB)[c->shift / LONG_BIT]);
		unsigned long lsb = htobe(((unsigned long *)c->CPB)[(c->shift + LONG_BIT - 1) / LONG_BIT]);
		buf = (msb << (c->shift % LONG_BIT)) | (lsb >> (-c->shift % LONG_BIT));
	}
	c->shift += v;
	c->ivlCurrRange <<= v;
	c->ivlOffset = (c->ivlOffset << v) | (buf >> (LONG_BIT - v));
}

static unsigned int get_ae(CABAC_ctx *c, uint8_t *state) {
	static const int rangeTabLps[4 * 64] = {
		128, 128, 128, 123, 116, 111, 105, 100, 95, 90, 85, 81, 77, 73, 69, 66,
		62, 59, 56, 53, 51, 48, 46, 43, 41, 39, 37, 35, 33, 32, 30, 29, 27, 26,
		24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 14, 13, 12, 12, 11, 11, 10,
		10, 9, 9, 8, 8, 7, 7, 7, 6, 6, 6, 2,
		176, 167, 158, 150, 142, 135, 128, 122, 116, 110, 104, 99, 94, 89, 85, 80,
		76, 72, 69, 65, 62, 59, 56, 53, 50, 48, 45, 43, 41, 39, 37, 35, 33, 31,
		30, 28, 27, 26, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 14, 13, 12,
		12, 11, 11, 10, 9, 9, 9, 8, 8, 7, 7, 2,
		208, 197, 187, 178, 169, 160, 152, 144, 137, 130, 123, 117, 111, 105, 100,
		95, 90, 86, 81, 77, 73, 69, 66, 63, 59, 56, 54, 51, 48, 46, 43, 41, 39,
		37, 35, 33, 32, 30, 29, 27, 26, 25, 23, 22, 21, 20, 19, 18, 17, 16, 15,
		15, 14, 13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 2,
		240, 227, 216, 205, 195, 185, 175, 166, 158, 150, 142, 135, 128, 122, 116,
		110, 104, 99, 94, 89, 85, 80, 76, 72, 69, 65, 62, 59, 56, 53, 50, 48, 45,
		43, 41, 39, 37, 35, 33, 31, 30, 28, 27, 25, 24, 23, 22, 21, 20, 19, 18,
		17, 16, 15, 14, 14, 13, 12, 12, 11, 11, 10, 9, 2,
	};
	static const int transIdx[2 * 128] = {
		0x7f, 0x7e, 0x4d, 0x4c, 0x4d, 0x4c, 0x4b, 0x4a, 0x4b, 0x4a, 0x4b, 0x4a,
		0x49, 0x48, 0x49, 0x48, 0x49, 0x48, 0x47, 0x46, 0x47, 0x46, 0x47, 0x46,
		0x45, 0x44, 0x45, 0x44, 0x43, 0x42, 0x43, 0x42, 0x43, 0x42, 0x41, 0x40,
		0x41, 0x40, 0x3f, 0x3e, 0x3d, 0x3c, 0x3d, 0x3c, 0x3d, 0x3c, 0x3b, 0x3a,
		0x3b, 0x3a, 0x39, 0x38, 0x37, 0x36, 0x37, 0x36, 0x35, 0x34, 0x35, 0x34,
		0x33, 0x32, 0x31, 0x30, 0x31, 0x30, 0x2f, 0x2e, 0x2d, 0x2c, 0x2d, 0x2c,
		0x2b, 0x2a, 0x2b, 0x2a, 0x27, 0x26, 0x27, 0x26, 0x25, 0x24, 0x25, 0x24,
		0x21, 0x20, 0x21, 0x20, 0x1f, 0x1e, 0x1f, 0x1e, 0x1b, 0x1a, 0x1b, 0x1a,
		0x19, 0x18, 0x17, 0x16, 0x17, 0x16, 0x13, 0x12, 0x13, 0x12, 0x11, 0x10,
		0x0f, 0x0e, 0x0d, 0x0c, 0x0b, 0x0a, 0x09, 0x08, 0x09, 0x08, 0x05, 0x04,
		0x05, 0x04, 0x03, 0x02, 0x01, 0x00, 0x00, 0x01,
		0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
		0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19,
		0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x22, 0x23, 0x24, 0x25,
		0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f, 0x30, 0x31,
		0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d,
		0x3e, 0x3f, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
		0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55,
		0x56, 0x57, 0x58, 0x59, 0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60, 0x61,
		0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d,
		0x6e, 0x6f, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79,
		0x7a, 0x7b, 0x7c, 0x7d, 0x7c, 0x7d, 0x7e, 0x7f,
	};
	
	unsigned int shift = LONG_BIT - 9 - __builtin_clzl(c->ivlCurrRange);
	unsigned int qRangeIdx = (c->ivlCurrRange >> c->shift) & 0xc0;
	unsigned long ivlLpsRange = (unsigned long)rangeTabLps[qRangeIdx | (*state >> 1)] << shift;
	unsigned long ivlMpsRange = c->ivlCurrRange - ivlLpsRange;
	unsigned long mask = ~(unsigned long)(c->ivlOffset - c->ivlCurrRange) >> (LONG_BIT - 1);
	c->ivlCurrRange = ivlMpsRange ^ ((ivlMpsRange ^ ivlLpsRange) & mask);
	c->ivlOffset -= ivlMpsRange & mask;
	int idx = *state ^ (int)mask;
	if (c->ivlCurrRange < 256)
		renorm(c, __builtin_clzl(c->ivlCurrRange) - 1);
	*state = transIdx[128 + idx];
	return idx & 1;
}

static inline unsigned int get_bypass(CABAC_ctx *c) {
	c->ivlCurrRange >>= 1;
	unsigned long mask = ~(unsigned long)(c->ivlOffset - c->ivlCurrRange) >> (LONG_BIT - 1);
	c->ivlOffset -= c->ivlCurrRange & mask;
	return -mask;
}



static const uint8_t ScanOrder4x4[3][16] = {
	{0, 4, 1, 8, 5, 2, 12, 9, 6, 3, 13, 10, 7, 14, 11, 15},
	{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
	{0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15},
};

static const uint8_t ScanOrder8x8[3][64] = {
	{0, 8, 1, 16, 9, 2, 24, 17, 10, 3, 32, 25, 18, 11, 4, 40, 33, 26, 19, 12, 5,
	48, 41, 34, 27, 20, 13, 6, 56, 49, 42, 35, 28, 21, 14, 7, 57, 50, 43, 36, 29,
	22, 15, 58, 51, 44, 37, 30, 23, 59, 52, 45, 38, 31, 60, 53, 46, 39, 61, 54,
	47, 62, 55, 63},
	{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
	21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
	40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
	59, 60, 61, 62, 63},
	{0, 8, 16, 24, 32, 40, 48, 56, 1, 9, 17, 25, 33, 41, 49, 57, 2, 10, 18, 26,
	34, 42, 50, 58, 3, 11, 19, 27, 35, 43, 51, 59, 4, 12, 20, 28, 36, 44, 52, 60,
	5, 13, 21, 29, 37, 45, 53, 61, 6, 14, 22, 30, 38, 46, 54, 62, 7, 15, 23, 31,
	39, 47, 55, 63},
};



struct Rage265_slice {
	Rage265_parameter_set p;
	unsigned int ctb_x:11;
	unsigned int ctb_y:11;
	unsigned int slice_type:2;
	unsigned int colour_plane_id:2;
	CABAC_ctx c;
};

#endif
