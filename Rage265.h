/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */
#ifndef RAGE265_H
#define RAGE265_H

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>

enum Rage265_errors {
	RAGE265_ERROR_NO_MEMORY          = 0,
	RAGE265_ERROR_PARSING_BITSTREAM  = 1,
	RAGE265_UNSUPPORTED_MULTIPLE_SPS = 2,
};

typedef struct Rage265_worker Rage265_worker;
typedef struct {
	unsigned int max_sub_layers:3;
	unsigned int temporal_id_nesting_flag:1;
	unsigned int ChromaArrayType:2;
	unsigned int separate_colour_plane_flag:1;
	unsigned int BitDepth_Y:4;
	unsigned int BitDepth_C:4;
	unsigned int log2_max_pic_order_cnt_lsb:5;
	unsigned int max_dec_pic_buffering:5;
	unsigned int max_num_reorder_pics:4;
	unsigned int MinCbLog2SizeY:3;
	unsigned int CtbLog2SizeY:3;
	unsigned int Log2MinTrafoSize:3;
	unsigned int Log2MaxTrafoSize:3;
	unsigned int max_transform_hierarchy_depth_inter:3;
	unsigned int max_transform_hierarchy_depth_intra:3;
	unsigned int amp_enabled_flag:1;
	unsigned int sample_adaptive_offset_enabled_flag:1;
	unsigned int PcmBitDepth_Y:4;
	unsigned int PcmBitDepth_C:4;
	unsigned int Log2MinIpcmCbSizeY:3;
	unsigned int Log2MaxIpcmCbSizeY:3;
	unsigned int pcm_loop_filter_disabled_flag:1;
	unsigned int temporal_mvp_enabled_flag:1;
	unsigned int strong_intra_smoothing_enabled_flag:1;
	uint16_t pic_width_in_luma_samples; // 15 significant bits
	uint16_t pic_height_in_luma_samples;
	uint16_t conf_win_left_offset; // in luma samples
	uint16_t conf_win_right_offset;
	uint16_t conf_win_top_offset;
	uint16_t conf_win_bottom_offset;
} Rage265_parameter_set;
typedef struct {
	uint8_t *image;
	unsigned int used_for_reference:1;
	unsigned int long_term_flag:1;
	int32_t PictureOrderCnt;
} Rage265_picture;
typedef struct {
	unsigned int max_workers;
	Rage265_worker *workers;
	pthread_mutex_t lock;
	pthread_cond_t worker_available;
	void *DPB;
	Rage265_parameter_set SPS;
	Rage265_parameter_set PPS[4];
} Rage265_ctx;

size_t Rage265_find_start_code(const uint8_t *buf, size_t len, unsigned int n);
unsigned int Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *buf, size_t len);

#endif
