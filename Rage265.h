/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */
#ifndef RAGE265_H
#define RAGE265_H

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>

enum Rage265_errors {
	RAGE265_ERROR_NO_MEMORY         = 0,
	RAGE265_ERROR_PARSING_BITSTREAM = 1,
};

typedef struct Rage265_worker Rage265_worker;
typedef struct {
	unsigned int vps_max_sub_layers:3;
	
} Rage265_VPS;
typedef struct {
	unsigned int max_workers;
	Rage265_worker *workers;
	pthread_mutex_t lock;
	pthread_cond_t worker_available;
	Rage265_VPS v;
} Rage265_ctx;

size_t Rage265_find_start_code(const uint8_t *buf, size_t len, unsigned int n);
unsigned int Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *buf, size_t len);

#endif
