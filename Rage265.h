/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */
#ifndef RAGE265_H
#define RAGE265_H

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>

enum Rage265_errors {
	RAGE265_ERROR_NO_MEMORY         = 0x1,
	RAGE265_ERROR_PARSING_BITSTREAM = 0x2,
};

typedef struct {
	unsigned int max_workers;
	void *workers;
	pthread_mutex_t lock;
	pthread_cond_t worker_available;
} Rage265_ctx;

size_t Rage265_find_start_code(const uint8_t *buf, size_t len, unsigned int n);
unsigned int Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *buf, size_t len);

#endif
