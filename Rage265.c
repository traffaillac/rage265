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
#include "Rage265.h"
#include "Rage265_common.h"

#include <stdlib.h>
#include <string.h>
#include <sys/sysinfo.h>



static void parse_scaling_list_data(Rage265_worker *w) {
	
}



static void parse_short_term_ref_pic_set(Rage265_worker *w) {
	
}



static void parse_vui_parameters(Rage265_worker *w) {
	
}



static void parse_profile_tier_level(Rage265_worker *w) {
	
}



static void parse_hrd_parameters(Rage265_worker *w, unsigned int cprms_present_flag, unsigned int max_sub_layers) {
	
}



/**
 * This function parses the SPS into a Rage265_parameter_set structure, and
 * stores it if no error was detected.
 */
static unsigned int parse_SPS(Rage265_ctx *r, Rage265_worker *w) {
	static const char * const chroma_format_idc_names[4] = {"4:0:0", "4:2:0", "4:2:2", "4:4:4"};
	
	Rage265_parameter_set s = {0};
	unsigned int sps_video_parameter_set_id = w->c.CPB[0] >> 4;
	s.max_sub_layers = min((w->c.CPB[0] >> 1) & 0x7, 6) + 1;
	s.temporal_id_nesting_flag = w->c.CPB[0] & 1;
	printf("<li>sps_video_parameter_set_id: <code>%u</code></li>\n"
		"<li>sps_max_sub_layers: <code>%u</code></li>\n"
		"<li>sps_temporal_id_nesting_flag: <code>%x</code></li>\n",
		sps_video_parameter_set_id,
		s.max_sub_layers,
		s.temporal_id_nesting_flag);
	w->c.shift = 8;
	parse_profile_tier_level(w);
	unsigned int sps_seq_parameter_set_id = get_ue(w->c.CPB, &w->c.shift, 15);
	unsigned int error_flags = (sps_seq_parameter_set_id > 0) << RAGE265_UNSUPPORTED_MULTIPLE_SPS;
	s.ChromaArrayType = get_ue(w->c.CPB, &w->c.shift, 3);
	printf("<li>sps_seq_parameter_set_id: <code>%u</code></li>\n"
		"<li>chroma_format_idc: <code>%u (%s)</code></li>\n",
		sps_seq_parameter_set_id,
		s.ChromaArrayType, chroma_format_idc_names[s.ChromaArrayType]);
	if (s.ChromaArrayType == 3) {
		s.separate_colour_plane_flag = get_u1(w->c.CPB, &w->c.shift);
		s.ChromaArrayType &= s.separate_colour_plane_flag - 1;
		printf("<li>separate_colour_plane_flag: <code>%x</code></li>\n",
			s.separate_colour_plane_flag);
	}
	s.pic_width_in_luma_samples = get_ue(w->c.CPB, &w->c.shift, 16888);
	s.pic_height_in_luma_samples = get_ue(w->c.CPB, &w->c.shift, 16888);
	unsigned int conformance_window_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>pic_width_in_luma_samples: <code>%u</code></li>\n"
		"<li>pic_height_in_luma_samples: <code>%u</code></li>\n",
		s.pic_width_in_luma_samples,
		s.pic_height_in_luma_samples);
	if (conformance_window_flag) {
		unsigned int shiftX = (s.ChromaArrayType == 1 || s.ChromaArrayType == 2);
		unsigned int shiftY = (s.ChromaArrayType == 1);
		s.conf_win_left_offset = min(get_ue16(w->c.CPB, &w->c.shift) << shiftX, s.pic_width_in_luma_samples);
		s.conf_win_right_offset = min(get_ue16(w->c.CPB, &w->c.shift) << shiftX, s.pic_width_in_luma_samples - s.conf_win_left_offset);
		s.conf_win_top_offset = min(get_ue16(w->c.CPB, &w->c.shift) << shiftY, s.pic_height_in_luma_samples);
		s.conf_win_bottom_offset = min(get_ue16(w->c.CPB, &w->c.shift) << shiftY, s.pic_height_in_luma_samples - s.conf_win_top_offset);
		printf("<li>conf_win_left_offset: <code>%u</code></li>\n"
			"<li>conf_win_right_offset: <code>%u</code></li>\n"
			"<li>conf_win_top_offset: <code>%u</code></li>\n"
			"<li>conf_win_bottom_offset: <code>%u</code></li>\n",
			s.conf_win_left_offset,
			s.conf_win_right_offset,
			s.conf_win_top_offset,
			s.conf_win_bottom_offset);
	}
	s.BitDepth_Y = get_ue(w->c.CPB, &w->c.shift, 6) + 8;
	s.BitDepth_C = get_ue(w->c.CPB, &w->c.shift, 6) + 8;
	s.log2_max_pic_order_cnt_lsb = get_ue(w->c.CPB, &w->c.shift, 12) + 4;
	unsigned int sps_sub_layer_ordering_info_present_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>BitDepth<sub>Y</sub>: <code>%u</code></li>\n"
		"<li>BitDepth<sub>C</sub>: <code>%u</code></li>\n"
		"<li>MaxPicOrderCntLsb: <code>%u</code></li>\n",
		s.BitDepth_Y,
		s.BitDepth_C,
		1 << s.log2_max_pic_order_cnt_lsb);
	for (unsigned int i = (s.max_sub_layers - 1) & -sps_sub_layer_ordering_info_present_flag; i < s.max_sub_layers; i++) {
		s.max_dec_pic_buffering = get_ue(w->c.CPB, &w->c.shift, 15) + 1;
		s.max_num_reorder_pics = get_ue(w->c.CPB, &w->c.shift, 15);
		unsigned int sps_max_latency_increase = get_ue(w->c.CPB, &w->c.shift, 4294967294) - 1;
		printf("<ul>\n"
			"<li>sps_max_dec_pic_buffering[%u]: <code>%u</code></li>\n"
			"<li>sps_max_num_reorder_pics[%u]: <code>%u</code></li>\n"
			"<li>sps_max_latency_increase[%u]: <code>%u</code></li>\n"
			"</ul>\n",
			i, s.max_dec_pic_buffering,
			i, s.max_num_reorder_pics,
			i, sps_max_latency_increase);
	}
	s.MinCbLog2SizeY = get_ue(w->c.CPB, &w->c.shift, 3) + 3;
	s.CtbLog2SizeY = min(max(s.MinCbLog2SizeY + get_ue8(w->c.CPB, &w->c.shift), 4), 6);
	s.Log2MinTrafoSize = min(get_ue8(w->c.CPB, &w->c.shift) + 2, min(s.MinCbLog2SizeY, 5));
	s.Log2MaxTrafoSize = min(s.Log2MinTrafoSize + get_ue8(w->c.CPB, &w->c.shift), min(s.CtbLog2SizeY, 5));
	s.max_transform_hierarchy_depth_inter = min(get_ue8(w->c.CPB, &w->c.shift), s.CtbLog2SizeY - s.Log2MinTrafoSize);
	s.max_transform_hierarchy_depth_intra = min(get_ue8(w->c.CPB, &w->c.shift), s.CtbLog2SizeY - s.Log2MinTrafoSize);
	unsigned int scaling_list_enabled_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>MinCbLog2SizeY: <code>%u</code></li>\n"
		"<li>CtbLog2SizeY: <code>%u</code></li>\n"
		"<li>Log2MinTrafoSize: <code>%u</code></li>\n"
		"<li>Log2MaxTrafoSize: <code>%u</code></li>\n"
		"<li>max_transform_hierarchy_depth_inter: <code>%u</code></li>\n"
		"<li>max_transform_hierarchy_depth_intra: <code>%u</code></li>\n"
		"<li>scaling_list_enabled_flag: <code>%x</code></li>\n",
		s.MinCbLog2SizeY,
		s.CtbLog2SizeY,
		s.Log2MinTrafoSize,
		s.Log2MaxTrafoSize,
		s.max_transform_hierarchy_depth_inter,
		s.max_transform_hierarchy_depth_intra,
		scaling_list_enabled_flag);
	if (scaling_list_enabled_flag) {
		unsigned int sps_scaling_list_data_present_flag = get_u1(w->c.CPB, &w->c.shift);
		if (sps_scaling_list_data_present_flag)
			parse_scaling_list_data(w);
	}
	s.amp_enabled_flag = get_u1(w->c.CPB, &w->c.shift);
	s.sample_adaptive_offset_enabled_flag = get_u1(w->c.CPB, &w->c.shift);
	unsigned int pcm_enabled_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>amp_enabled_flag: <code>%x</code></li>\n"
		"<li>sample_adaptive_offset_enabled_flag: <code>%x</code></li>\n"
		"<li>pcm_enabled_flag: <code>%x</code></li>\n",
		s.amp_enabled_flag,
		s.sample_adaptive_offset_enabled_flag,
		pcm_enabled_flag);
	if (pcm_enabled_flag) {
		s.PcmBitDepth_Y = min(get_uv(w->c.CPB, &w->c.shift, 4) + 1, s.BitDepth_Y);
		s.PcmBitDepth_C = min(get_uv(w->c.CPB, &w->c.shift, 4) + 1, s.BitDepth_C);
		s.Log2MinIpcmCbSizeY = min(max(get_ue8(w->c.CPB, &w->c.shift) + 3, s.MinCbLog2SizeY), min(s.CtbLog2SizeY, 5));
		s.Log2MaxIpcmCbSizeY = min(s.Log2MinIpcmCbSizeY + get_ue8(w->c.CPB, &w->c.shift), min(s.CtbLog2SizeY, 5));
		s.pcm_loop_filter_disabled_flag = get_u1(w->c.CPB, &w->c.shift);
		printf("<li>PcmBitDepth<sub>Y</sub>: <code>%u</code></li>\n"
			"<li>PcmBitDepth<sub>C</sub>: <code>%u</code></li>\n"
			"<li>Log2MinIpcmCbSizeY: <code>%u</code></li>\n"
			"<li>Log2MaxIpcmCbSizeY: <code>%u</code></li>\n"
			"<li>pcm_loop_filter_disabled_flag: <code>%x</code></li>\n",
			s.PcmBitDepth_Y,
			s.PcmBitDepth_C,
			s.Log2MinIpcmCbSizeY,
			s.Log2MaxIpcmCbSizeY,
			s.pcm_loop_filter_disabled_flag);
	}
	unsigned int num_short_term_ref_pic_sets = get_ue(w->c.CPB, &w->c.shift, 64);
	printf("<li>num_short_term_ref_pic_sets: <code>%u</code></li>\n",
		num_short_term_ref_pic_sets);
	for (unsigned int i = 0; i < num_short_term_ref_pic_sets; i++)
		parse_short_term_ref_pic_set(w);
	unsigned int long_term_ref_pics_present_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>long_term_ref_pics_present_flag: <code>%x</code></li>\n",
		long_term_ref_pics_present_flag);
	if (long_term_ref_pics_present_flag) {
		unsigned int num_long_term_ref_pics_sps = get_ue(w->c.CPB, &w->c.shift, 32);
		printf("<li>num_long_term_ref_pics_sps: <code>%u</code></li>\n",
			num_long_term_ref_pics_sps);
		for (unsigned int i = 0; i < num_long_term_ref_pics_sps; i++) {
			unsigned int lt_ref_pic_poc_lsb_sps = get_uv(w->c.CPB, &w->c.shift, s.log2_max_pic_order_cnt_lsb);
			unsigned int used_by_curr_pic_lt_sps_flag = get_u1(w->c.CPB, &w->c.shift);
			printf("<ul>\n"
				"<li>lt_ref_pic_poc_lsb_sps[%u]: <code>%u</code></li>\n"
				"<li>used_by_curr_pic_lt_sps_flag[%u]: <code>%u</code></li>\n"
				"</ul>\n",
				i, lt_ref_pic_poc_lsb_sps,
				i, used_by_curr_pic_lt_sps_flag);
		}
	}
	s.temporal_mvp_enabled_flag = get_u1(w->c.CPB, &w->c.shift);
	s.strong_intra_smoothing_enabled_flag = get_u1(w->c.CPB, &w->c.shift);
	unsigned int vui_parameters_present_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>sps_temporal_mvp_enabled_flag: <code>%x</code></li>\n"
		"<li>strong_intra_smoothing_enabled_flag: <code>%x</code></li>\n",
		s.temporal_mvp_enabled_flag,
		s.strong_intra_smoothing_enabled_flag);
	if (vui_parameters_present_flag)
		parse_vui_parameters(w);
	unsigned int sps_extension_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>sps_extension_flag: <code>%x</code></li>\n", sps_extension_flag);
	if (sps_extension_flag && w->c.shift < w->c.lim)
		w->c.shift = w->c.lim;
	error_flags |= (w->c.shift != w->c.lim) << RAGE265_ERROR_PARSING_BITSTREAM;
	if (error_flags)
		return error_flags;
	
	/* Clear the CPBs and reallocate the DPB when the image format changes. */
	if ((s.ChromaArrayType ^ r->SPS.ChromaArrayType) |
		(s.pic_width_in_luma_samples ^ r->SPS.pic_width_in_luma_samples) |
		(s.pic_height_in_luma_samples ^ r->SPS.pic_height_in_luma_samples) |
		(s.BitDepth_Y ^ r->SPS.BitDepth_Y) | (s.BitDepth_C ^ r->SPS.BitDepth_C) |
		(s.max_dec_pic_buffering ^ r->SPS.max_dec_pic_buffering)) {
		pthread_mutex_lock(&r->lock);
		for (unsigned int i = 0; i < r->max_workers; i++) {
			while (r->workers[i].target != NULL)
				pthread_cond_wait(&r->workers[i].target_changed, &r->lock);
			if (r->workers[i].c.CPB != NULL) {
				free((uint8_t *)r->workers[i].c.CPB);
				r->workers[i].c.CPB = NULL;
				r->workers[i].CPB_size = 0;
			}
		}
		pthread_mutex_unlock(&r->lock);
		s.pic_width_in_luma_samples &= -1 << s.MinCbLog2SizeY;
		s.pic_height_in_luma_samples &= -1 << s.MinCbLog2SizeY;
		unsigned int PicSizeInSamplesY = s.pic_width_in_luma_samples * s.pic_height_in_luma_samples;
		if (PicSizeInSamplesY > 35651584) {
			s.pic_height_in_luma_samples = 35651584 / s.pic_width_in_luma_samples;
			PicSizeInSamplesY = s.pic_width_in_luma_samples * s.pic_height_in_luma_samples;
		}
		unsigned int PicSizeInSamplesC = PicSizeInSamplesY / 4 * (1 << s.ChromaArrayType >> 1);
		size_t luma_size = PicSizeInSamplesY << ((s.BitDepth_Y - 1) / 8);
		size_t chroma_size = PicSizeInSamplesC << ((s.BitDepth_C - 1) / 8);
		size_t picture_size = sizeof(Rage265_picture) + luma_size + 2 * chroma_size;
		if (r->DPB != NULL)
			free(r->DPB);
		r->DPB = malloc(picture_size * (s.max_dec_pic_buffering - 1 + r->max_workers));
		for (unsigned int i = 0; i < s.max_dec_pic_buffering - 1 + r->max_workers; i++) {
			Rage265_picture *p = r->DPB + i * picture_size;
			p->image = (uint8_t *)p + sizeof(*p);
			p->used_for_reference = 0;
			p->long_term_flag = 0;
			p->PictureOrderCnt = INT32_MIN;
		}
	}
	r->SPS = s;
	return 0;
}



/**
 * This function currently only prints the VPS to stdout.
 */
static unsigned int parse_VPS(Rage265_ctx *r, Rage265_worker *w) {
	unsigned int u = htobe16(*(uint16_t *)w->c.CPB);
	unsigned int vps_video_parameter_set_id = u >> 12;
	unsigned int vps_max_layers = ((u >> 4) & 0x3f) + 1;
	unsigned int vps_max_sub_layers = min((u >> 1) & 0x7, 6) + 1;
	unsigned int vps_temporal_id_nesting_flag = u & 1;
	printf("<li>vps_video_parameter_set_id: <code>%u</code></li>\n"
		"<li>vps_max_layers: <code>%u</code></li>\n"
		"<li>vps_max_sub_layers: <code>%u</code></li>\n"
		"<li>vps_temporal_id_nesting_flag: <code>%x</code></li>\n",
		vps_video_parameter_set_id,
		vps_max_layers,
		vps_max_sub_layers,
		vps_temporal_id_nesting_flag);
	w->c.shift = 32;
	parse_profile_tier_level(w);
	unsigned int vps_sub_layer_ordering_info_present_flag = get_u1(w->c.CPB, &w->c.shift);
	for (unsigned int i = (vps_max_sub_layers - 1) & -vps_sub_layer_ordering_info_present_flag; i < vps_max_sub_layers; i++) {
		unsigned int vps_max_dec_pic_buffering = get_ue(w->c.CPB, &w->c.shift, 15) + 1;
		unsigned int vps_max_num_reorder_pics = get_ue(w->c.CPB, &w->c.shift, 15);
		unsigned int vps_max_latency_increase = get_ue(w->c.CPB, &w->c.shift, 4294967294) - 1;
		printf("<ul>\n"
			"<li>vps_max_dec_pic_buffering[%u]: <code>%u</code></li>\n"
			"<li>vps_max_num_reorder_pics[%u]: <code>%u</code></li>\n"
			"<li>vps_max_latency_increase[%u]: <code>%u</code></li>\n"
			"</ul>\n",
			i, vps_max_dec_pic_buffering,
			i, vps_max_num_reorder_pics,
			i, vps_max_latency_increase);
	}
	unsigned int vps_max_layer_id = get_uv(w->c.CPB, &w->c.shift, 6);
	unsigned int vps_num_layer_sets = get_ue(w->c.CPB, &w->c.shift, 1023) + 1;
	printf("<li>vps_max_layer_id: <code>%u</code></li>\n"
		"<li>vps_num_layer_sets: <code>%u</code></li>\n",
		vps_max_layer_id,
		vps_num_layer_sets);
	for (unsigned int i = 1; i < vps_num_layer_sets; i++) {
		for (unsigned int j = 0; j <= vps_max_layer_id; j++) {
			unsigned int layer_id_included_flag = get_u1(w->c.CPB, &w->c.shift);
			printf("<ul><li>layer_id_included_flag[%u][%u]: <code>%x</code></li></ul>\n",
				i, j, layer_id_included_flag);
		}
	}
	unsigned int vps_timing_info_present_flag = get_u1(w->c.CPB, &w->c.shift);
	if (vps_timing_info_present_flag) {
		unsigned int vps_num_units_in_tick = get_uv(w->c.CPB, &w->c.shift, 32);
		unsigned int vps_time_scale = get_uv(w->c.CPB, &w->c.shift, 32);
		unsigned int vps_poc_proportional_to_timing_flag = get_u1(w->c.CPB, &w->c.shift);
		printf("<li>vps_num_units_in_tick: <code>%u</code></li>\n"
			"<li>vps_time_scale: <code>%u</code></li>\n"
			"<li>vps_poc_proportional_to_timing_flag: <code>%x</code></li>\n",
			vps_num_units_in_tick,
			vps_time_scale,
			vps_poc_proportional_to_timing_flag);
		if (vps_poc_proportional_to_timing_flag) {
			unsigned int vps_num_ticks_poc_diff_one = get_ue(w->c.CPB, &w->c.shift, 4294967294) + 1;
			printf("<li>vps_num_ticks_poc_diff_one: <code>%u</code></li>\n",
				vps_num_ticks_poc_diff_one);
		}
		unsigned int vps_num_hrd_parameters = get_ue(w->c.CPB, &w->c.shift, 1024);
		printf("<li>vps_num_hrd_parameters: <code>%u</code></li>\n",
			vps_num_hrd_parameters);
		for (unsigned int i = 0; i < vps_num_hrd_parameters; i++) {
			unsigned int hrd_layer_set_idx = get_ue(w->c.CPB, &w->c.shift, 1023);
			printf("<ul>\n"
				"<li>hrd_layer_set_idx[%u]: <code>%u</code></li>\n",
				i, hrd_layer_set_idx);
			unsigned int cprms_present_flag = 0;
			if (i > 0)
				cprms_present_flag = get_u1(w->c.CPB, &w->c.shift);
			parse_hrd_parameters(w, cprms_present_flag, vps_max_sub_layers);
			printf("</ul>\n");
		}
	}
	unsigned int vps_extension_flag = get_u1(w->c.CPB, &w->c.shift);
	printf("<li>vps_extension_flag: <code>%x</code></li>\n", vps_extension_flag);
	if (vps_extension_flag && w->c.shift < w->c.lim)
		w->c.shift = w->c.lim;
	return (w->c.shift != w->c.lim) << RAGE265_ERROR_PARSING_BITSTREAM;
}



/**
 * Find the start of the next 00 00 0n pattern, returning len if none was found.
 */
#ifdef __SSSE3__
size_t Rage265_find_start_code(const uint8_t *buf, size_t len, unsigned int n) {
	ptrdiff_t chunk = (uint8_t *)((uintptr_t)buf & -sizeof(__m128i)) - buf;
	for (size_t u = 0; chunk < (ptrdiff_t)len; u = chunk += sizeof(__m128i)) {
		/* Skip chunks without a zero odd byte. */
		if (__builtin_expect((_mm_movemask_epi8(_mm_cmpeq_epi8(*(__m128i *)(buf + chunk), _mm_setzero_si128())) & 0xaaaa) != 0, 0)) {
			size_t lim = min(chunk + sizeof(__m128i) + 2, len);
			for (unsigned int start_code = -1; u < lim; u++) {
				start_code = ((start_code & 0xffff) << 8) | buf[u];
				if (start_code == n)
					return u - 2;
			}
		}
	}
	return len;
}
#endif



/**
 * Parse a NAL unit, returning a bitfield of error codes.
 */
unsigned int Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *buf, size_t len) {
	static const char * const nal_unit_type_names[64] = {
		[0] = "Trailing non-reference picture",
		[1] = "Trailing reference picture",
		[2] = "Temporal Sub-layer Access non-reference picture",
		[3] = "Temporal Sub-layer Access reference picture",
		[4] = "Step-wise Temporal Sub-layer Access non-reference picture",
		[5] = "Step-wise Temporal Sub-layer Access reference picture",
		[6] = "Random Access Decodable Leading non-reference picture",
		[7] = "Random Access Decodable Leading reference picture",
		[8] = "Random Access Skipped Leading non-reference picture",
		[9] = "Random Access Skipped Leading reference picture",
		[10 ... 15] = "unknown",
		[16] = "Broken Link Access picture with Leading Picture",
		[17] = "Broken Link Access picture with Random Access Decodable Leading picture",
		[18] = "Broken Link Access picture without Leading Picture",
		[19] = "Instantaneous Decoding Refresh picture with Random Access Decodable Leading picture",
		[20] = "Instantaneous Decoding Refresh picture without Leading Picture",
		[21] = "Clean Random Access picture",
		[22 ... 31] = "unknown",
		[32] = "Video Parameter Set",
		[33] = "Sequence Parameter Set",
		[34] = "Picture Parameter Set",
		[35] = "Access Unit Delimiter",
		[36] = "End Of Sequence",
		[37] = "End Of Bitstream",
		[38] = "Filler Data",
		[39] = "Prefix Supplemental Enhancement Information",
		[40] = "Suffix Supplemental Enhancement Information",
		[41 ... 63] = "unknown",
	};
	typedef unsigned int (*Parser)(Rage265_ctx *, Rage265_worker *);
	static const Parser parse_nal_unit[64] = {
		[32] = parse_VPS,
		[33] = parse_SPS,
	};
	
	/* On first call, initialise the main structure. */
	if (r->max_workers == 0)
		r->max_workers = get_nprocs();
	if (r->workers == NULL) {
		r->workers = calloc(r->max_workers, sizeof(*r->workers));
		if (r->workers == NULL)
			return 1 << RAGE265_ERROR_NO_MEMORY;
		r->lock = (pthread_mutex_t)PTHREAD_MUTEX_INITIALIZER;
		r->worker_available = (pthread_cond_t)PTHREAD_COND_INITIALIZER;
	}
	
	/* Assign an available worker to this NAL payload. */
	pthread_mutex_lock(&r->lock);
	unsigned int i = r->max_workers;
	while (i == r->max_workers) {
		for (i = 0; i < r->max_workers && r->workers[i].target != NULL; i++);
		if (i == r->max_workers)
			pthread_cond_wait(&r->worker_available, &r->lock);
	}
	Rage265_worker *w = r->workers + i;
	pthread_mutex_unlock(&r->lock);
	
	/* Allocate the CPB to let the worker process it asynchronously. */
	if (len < 2)
		return 1 << RAGE265_ERROR_PARSING_BITSTREAM;
	const unsigned int suffix_size = 128;
	size_t CPB_size = len - 2 + suffix_size;
	if (CPB_size > 800000000 / 8 + suffix_size) { // Level 6.2, High tier
		CPB_size = 800000000 / 8 + suffix_size;
		len = 800000000 / 8;
	}
	if (w->CPB_size < CPB_size) {
		w->CPB_size = CPB_size;
		if (w->c.CPB != NULL)
			free((uint8_t *)w->c.CPB);
		w->c.CPB = malloc(CPB_size);
		if (w->c.CPB == NULL)
			return 1 << RAGE265_ERROR_NO_MEMORY;
	}
	
	/* Copy the CPB while removing every emulation_prevention_three_byte. */
	uint8_t *dst = (uint8_t *)w->c.CPB;
	unsigned int u = 2;
	while (u < len) {
		unsigned int copy = Rage265_find_start_code(buf + u, len - u, 3);
		memcpy(dst, buf + u, copy + 2 * (copy < len - u));
		dst += copy + 2;
		u += copy + 3;
	}
	dst -= 3; // Skips one cabac_zero_word as a side-effect
	
	/* Trim every cabac_zero_word, delimit the SODB, and append the safety suffix. */
	while (*dst == 0)
		dst -= 2;
	w->c.shift = 0;
	w->c.lim = 8 * (dst - w->c.CPB) + 7 - __builtin_ctz(*dst);
	memset(dst + 1, 0xff, suffix_size);
	
	/* Parse the nal_unit_header(). */
	unsigned int nal_unit_header = (buf[0] << 8) | buf[1];
	w->nal_unit_type = nal_unit_header >> 9;
	unsigned int nuh_layer_id = (nal_unit_header >> 3) & 0x3f;
	unsigned int nuh_temporal_id = (nal_unit_header - 1) & 0x7;
	printf("<ul class=\"frame\">\n"
		"<li%s>nal_unit_type: <code>%u (%s)</code></li>\n"
		"<li>nuh_layer_id: <code>%u</code></li>\n"
		"<li>nuh_temporal_id: <code>%u</code></li>\n",
		red_if(parse_nal_unit[w->nal_unit_type] == NULL), w->nal_unit_type, nal_unit_type_names[w->nal_unit_type],
		nuh_layer_id,
		nuh_temporal_id);
	
	/* Branch on nal_unit_type. */
	unsigned int error_flags = 0;
	if (nuh_layer_id == 0 && parse_nal_unit[w->nal_unit_type] != NULL)
		error_flags = parse_nal_unit[w->nal_unit_type](r, w);
	if (error_flags)
		printf("<li style=\"color: red\">Error 0x%x</li>\n", error_flags);
	printf("</ul>\n");
	return error_flags;
}
