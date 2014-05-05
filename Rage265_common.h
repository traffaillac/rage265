/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */
#ifndef RAGE265_COMMON_H
#define RAGE265_COMMON_H

#include <assert.h>
#include <endian.h>
#include <limits.h>
#include <stdint.h>

#ifndef WORD_BIT
#define WORD_BIT (sizeof(int) * 8)
#endif
#ifndef LONG_BIT
#define LONG_BIT (sizeof(long) * 8)
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

static inline __attribute__((always_inline)) unsigned int get_ue(const uint8_t *CPB, unsigned int *shift, unsigned int upper) {
	assert(upper<4294967296);
	unsigned int leadingZeroBits, res;
	if (upper <= 31) {
		uint16_t buf = (((CPB[*shift / 8] << 8) | CPB[*shift / 8 + 1]) << (*shift % 8));
		leadingZeroBits = __builtin_clz(buf | 0x0400) - WORD_BIT + 16;
		res = (buf >> (16 - (2 * leadingZeroBits + 1))) - 1;
	} else if (upper <= 65534) {
		unsigned int msb = be32toh(((uint32_t *)CPB)[*shift / 32]);
		unsigned int lsb = be32toh(((uint32_t *)CPB)[(*shift + 31) / 32]);
		uint32_t buf = (msb << (*shift % 32)) | (lsb >> (-*shift % 32));
		leadingZeroBits = __builtin_clz(buf | 0x00010000) - WORD_BIT + 32;
		res = (buf >> (32 - (2 * leadingZeroBits + 1))) - 1;
	} else {
		uint64_t msb = be64toh(((uint64_t *)CPB)[*shift / 64]);
		uint64_t lsb = be64toh(((uint64_t *)CPB)[(*shift + 63) / 64]);
		uint64_t buf = (msb << (*shift % 64)) | (lsb >> (-*shift % 64));
		leadingZeroBits = clz64(buf | 0x0000000100000000);
		res = (buf >> (64 - (2 * leadingZeroBits + 1))) - 1;
	}
	*shift += 2 * leadingZeroBits + 1;
	return (res < upper) ? res : upper;
}

static inline __attribute__((always_inline)) int get_se(const uint8_t *CPB, unsigned int *shift, int lower, int upper) {
	unsigned int codeNum = get_ue(CPB, shift, max(-lower * 2, upper * 2 - 1));
	unsigned int abs = (codeNum+ 1) / 2;
	unsigned int sign = (codeNum % 2) - 1;
	int res = (abs ^ sign) - sign;
	return (-lower * 2 < upper * 2 - 1) ? max(res, lower) : (-lower * 2 > upper * 2 - 1) ? min(res, upper) : res;
}

static inline __attribute__((always_inline)) unsigned int get_uv(const uint8_t *CPB, unsigned int *shift, unsigned int v) {
	unsigned int msb = be32toh(((uint32_t *)CPB)[*shift / 32]);
	unsigned int lsb = be32toh(((uint32_t *)CPB)[(*shift + 31) / 32]);
	uint32_t buf = (msb << (*shift % 32)) | (lsb >> (-*shift % 32));
	*shift += v;
	return buf >> (32 - v);
}

static inline __attribute__((always_inline)) unsigned int get_u1(const uint8_t *CPB, unsigned int *shift) {
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
#define betoh be32toh
#elif ULONG_MAX == 18446744073709551615U
#define betoh be64toh
#endif

static inline void renorm(CABAC_ctx *c, unsigned int v) {
	unsigned long buf = 0;
	if (c->shift < c->lim) {
		unsigned long msb = betoh(((unsigned long *)c->CPB)[c->shift / LONG_BIT]);
		unsigned long lsb = betoh(((unsigned long *)c->CPB)[(c->shift + LONG_BIT - 1) / LONG_BIT]);
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



typedef struct {
	pthread_t thread_id;
	pthread_mutex_t *lock;
	pthread_cond_t target_changed;
	pthread_cond_t *target_complete;
	uint8_t *target;
	CABAC_ctx c;
	unsigned int CPB_size; // in bytes, 27 significant bits
} Worker_ctx;

#endif
