/**
 * Copyright (c) 2014 Thibault Raffaillac <traf@kth.se>
 */
#ifndef RAGE265_H
#define RAGE265_H

#include <stddef.h>
#include <stdint.h>

enum Rage265_status {
	RAGE265_ERROR_NO_MEMORY,
	RAGE265_UNSUPPORTED_MULTIPLE_SPS,
	RAGE265_UNSUPPORTED_MORE_THAN_FOUR_PPS,
	RAGE265_UNSUPPORTED_DEPENDENT_SLICES,
};

typedef struct Rage265_worker Rage265_worker;
typedef struct {
	unsigned int max_sub_layers:3;
	unsigned int temporal_id_nesting_flag:1;
	unsigned int general_profile_space:2;
	unsigned int general_progressive_source_flag:1;
	unsigned int general_interlaced_source_flag:1;
	unsigned int ChromaArrayType:2;
	unsigned int separate_colour_plane_flag:1;
	unsigned int BitDepth_Y:4;
	unsigned int BitDepth_C:4;
	unsigned int QpBdOffset_Y:6;
	unsigned int QpBdOffset_C:6;
	unsigned int log2_max_pic_order_cnt_lsb:5;
	unsigned int max_dec_pic_buffering:5;
	unsigned int max_num_reorder_pics:4;
	unsigned int MinCbLog2SizeY:3;
	unsigned int CtbLog2SizeY:3;
	unsigned int PicWidthInCtbsY:11;
	unsigned int PicHeightInCtbsY:11;
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
	unsigned int num_short_term_ref_pic_sets:7;
	unsigned int long_term_ref_pics_present_flag:1;
	unsigned int num_long_term_ref_pics_sps:6;
	unsigned int temporal_mvp_enabled_flag:1;
	unsigned int strong_intra_smoothing_enabled_flag:1;
	unsigned int dependent_slice_segments_enabled_flag:1;
	unsigned int output_flag_present_flag:1;
	unsigned int num_extra_slice_header_bits:3;
	unsigned int sign_data_hiding_enabled_flag:1;
	unsigned int cabac_init_present_flag:1;
	unsigned int constrained_intra_pred_flag:1;
	unsigned int transform_skip_enabled_flag:1;
	unsigned int cu_qp_delta_enabled_flag:1;
	unsigned int Log2MinCuQpDeltaSize:3;
	int cb_qp_offset:5;
	int cr_qp_offset:5;
	unsigned int pps_slice_chroma_qp_offsets_present_flag:1;
	unsigned int weighted_pred_flags:2; // weighted_pred_flag << 1 | weighted_bipred_flag
	unsigned int transquant_bypass_enabled_flag:1;
	unsigned int entropy_coding_sync_enabled_flag:1;
	unsigned int num_tile_columns:5;
	unsigned int num_tile_rows:5;
	unsigned int loop_filter_across_tiles_enabled_flag:1;
	unsigned int loop_filter_across_slices_enabled_flag:1;
	unsigned int deblocking_filter_override_enabled_flag:1;
	unsigned int deblocking_filter_disabled_flag:1;
	int beta_offset:4;
	int tc_offset:4;
	unsigned int lists_modification_present_flag:1;
	unsigned int Log2ParMrgLevel:3;
	unsigned int slice_segment_header_extension_present_flag:1;
	uint16_t pic_width_in_luma_samples; // 15 significant bits
	uint16_t pic_height_in_luma_samples;
	uint16_t conf_win_left_offset; // in luma samples
	uint16_t conf_win_right_offset;
	uint16_t conf_win_top_offset;
	uint16_t conf_win_bottom_offset;
	uint16_t lt_ref_pic_poc_lsb_sps[32];
	uint32_t used_by_curr_pic_lt_sps_flags;
	uint16_t colBd[23]; // 11 significant bits
	uint16_t rowBd[21];
	uint8_t num_ref_idx_active[2]; // 4 significant bits each
	int8_t Qp; // 7 significant bits
	uint8_t ScalingFactor4x4[6][16] __attribute__((aligned));
	uint8_t ScalingFactor8x8[6][64] __attribute__((aligned));
	uint8_t ScalingFactor16x16[6][64] __attribute__((aligned));
	uint8_t ScalingFactor32x32[2][64] __attribute__((aligned));
} Rage265_parameter_set;
typedef struct {
	uint8_t *image;
	unsigned int used_for_reference:1;
	unsigned int long_term_flag:1;
	int32_t PicOrderCntVal;
} Rage265_picture;
typedef struct {
	uint8_t *CPB;
	unsigned int CPB_size; // in bytes, 27 significant bits
	unsigned int nal_unit_type:6;
	int32_t prevPicOrderCntVal;
	Rage265_picture *DPB;
	Rage265_parameter_set SPS;
	Rage265_parameter_set PPSs[4];
	uint16_t short_term_RPSs[64][16]; // [15] = NumDeltaPocs << 8 | NumNegativePics,
		// [0..14] = abs_delta_poc_minus1 << 1 | used_by_curr_pic_flag,
		// sorted by ascending values of DeltaPoc (unlike spec!).
} Rage265_ctx;

size_t Rage265_find_start_code(const uint8_t *buf, size_t len, unsigned int n);
const Rage265_picture *Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *buf, size_t len);

#endif
