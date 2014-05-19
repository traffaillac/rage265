// TODO: Deblocking is carried in 3 dedicated threads with row callbacks executing before the image is released for output
// Try ORing instead of setting every parameter_set field to see if generated code is smaller

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

#include <stdlib.h>
#include <string.h>



static const uint8_t intra_ScalingFactor[64] __attribute__((aligned)) = {
	16, 16, 16, 16, 17, 18, 21, 24,
	16, 16, 16, 16, 17, 19, 22, 25,
	16, 16, 17, 18, 20, 22, 25, 29,
	16, 16, 18, 21, 24, 27, 31, 36,
	17, 17, 20, 24, 30, 35, 41, 47,
	18, 19, 22, 27, 35, 44, 54, 65,
	21, 22, 25, 31, 41, 54, 70, 88,
	24, 25, 29, 36, 47, 65, 88, 115,
};

static const uint8_t inter_ScalingFactor[64] __attribute__((aligned)) = {
	16, 16, 16, 16, 17, 18, 20, 24,
	16, 16, 16, 17, 18, 20, 24, 25,
	16, 16, 17, 18, 20, 24, 25, 28,
	16, 17, 18, 20, 24, 25, 28, 33,
	17, 18, 20, 24, 25, 28, 33, 41,
	18, 20, 24, 25, 28, 33, 41, 54,
	20, 24, 25, 28, 33, 41, 54, 71,
	24, 25, 28, 33, 41, 54, 71, 91,
};



/**
 * This function parses a single short_term_ref_pic_set(), fills the proper
 * st_RPS entry, and returns the last shift value.
 */
static unsigned int parse_short_term_ref_pic_set(uint16_t (*st_RPS)[16],
	unsigned int stRpsIdx, const uint8_t *CPB, unsigned int shift,
	unsigned int num_short_term_ref_pic_sets, const Rage265_parameter_set *p)
{
	unsigned int inter_ref_pic_set_prediction_flag = 0;
	if (stRpsIdx > 0)
		inter_ref_pic_set_prediction_flag = get_u1(CPB, &shift);
	if (inter_ref_pic_set_prediction_flag) {
		unsigned int RefRpsIdx = stRpsIdx - 1;
		if (stRpsIdx == num_short_term_ref_pic_sets)
			RefRpsIdx -= min(get_ue32(CPB, &shift), RefRpsIdx);
		unsigned int NumNegativePics = st_RPS[RefRpsIdx][15] & 0xff;
		unsigned int NumDeltaPocs = st_RPS[RefRpsIdx][15] >> 8;
		unsigned int delta_rps_sign = get_u1(CPB, &shift);
		unsigned int abs_delta_rps = get_ue(CPB, &shift, 32767) + 1;
		int deltaRps = (abs_delta_rps ^ -delta_rps_sign) + delta_rps_sign;
		unsigned int used_by_curr_pic_flags = 0;
		unsigned int use_delta_flags = 0x1fffe;
		for (unsigned int j = 0; j <= NumDeltaPocs; j++) {
			unsigned int reorder = (j < NumNegativePics) ? NumNegativePics - j :
				(j < NumDeltaPocs) ? j + 2 : (NumNegativePics + 1);
			unsigned int used_by_curr_pic_flag = get_u1(CPB, &shift);
			used_by_curr_pic_flags |= used_by_curr_pic_flag << reorder;
			if (!used_by_curr_pic_flag)
				use_delta_flags ^= (get_u1(CPB, &shift) ^ 1) << reorder;
		}
		unsigned int i = 0, k = 0;
		for (unsigned int j = 0; j <= NumDeltaPocs; j++) {
			int dPoc = deltaRps;
			if (j < NumNegativePics)
				dPoc -= (st_RPS[RefRpsIdx][j] >> 1) + 1;
			else if (j > NumNegativePics)
				dPoc += (st_RPS[RefRpsIdx][j - 1] >> 1) + 1;
			if (((use_delta_flags >>= 1) & 1) && dPoc != 0) {
				st_RPS[stRpsIdx][i++] = (min(abs(dPoc) - 1, 32767) << 1) |
					((used_by_curr_pic_flags >>= 1) & 1);
				if (dPoc < 0)
					k = i;
			}
		}
		st_RPS[stRpsIdx][15] = (i << 8) | k;
	} else {
		unsigned int num_negative_pics = min(get_ue8(CPB, &shift),
			p->max_dec_pic_buffering - 1);
		unsigned int NumDeltaPocs = min(num_negative_pics + get_ue8(CPB, &shift),
			p->max_dec_pic_buffering - 1);
		st_RPS[stRpsIdx][15] = (NumDeltaPocs << 8) | num_negative_pics;
		int delta_poc_minus1 = -1;
		for (unsigned int i = 0; i < NumDeltaPocs; i++) {
			if (i == num_negative_pics)
				delta_poc_minus1 = -1;
			delta_poc_minus1 = min(delta_poc_minus1 + get_ue32(CPB, &shift) + 1, 32767);
			unsigned int used_by_curr_pic_flag = get_u1(CPB, &shift);
			unsigned int j = (i < num_negative_pics) ? num_negative_pics - i - 1 : i;
			st_RPS[stRpsIdx][j] = (delta_poc_minus1 << 1) | used_by_curr_pic_flag;
		}
	}
	printf("<ul>\n");
	for (unsigned int i = 0; i < st_RPS[stRpsIdx][15] >> 8; i++) {
		unsigned int abs_delta_poc = (st_RPS[stRpsIdx][i] >> 1) + 1;
		printf("<li>DeltaPoc[%u][%u]: <code>%d, %s</code></li>\n",
			stRpsIdx, i, (i < (st_RPS[stRpsIdx][15] & 0xff)) ? -abs_delta_poc : abs_delta_poc, (st_RPS[stRpsIdx][i] & 1) ? "used" : "follow");
	}
	printf("</ul>\n");
	return shift;
}



static inline unsigned int parse_ref_pic_list_modification(const uint8_t *CPB,
	unsigned int shift)
{
	
	return shift;
}



static inline unsigned int parse_pred_weight_table(const uint8_t *CPB,
	unsigned int shift)
{
	
	return shift;
}



/**
 * This function parses a section_segment_header(), and calls
 * parse_slice_segment_data() if no error was detected.
 */
static void parse_slice_segment_header(Rage265_ctx *r, unsigned int lim) {
	static const char * const slice_type_names[3] = {"B", "P", "I"};
	static const char * const colour_plane_id_names[3] = {"Y", "Cb", "Cr"};
	
	unsigned int shift = 0;
	unsigned int first_slice_segment_in_pic_flag = get_u1(r->CPB, &shift);
	printf("<li>first_slice_segment_in_pic_flag: <code>%x</code></li>\n",
		first_slice_segment_in_pic_flag);
	if (r->nal_unit_type >= 16 && r->nal_unit_type <= 23) {
		unsigned int no_output_of_prior_pics_flag = get_u1(r->CPB, &shift);
		printf("<li>no_output_of_prior_pics_flag: <code>%x</code></li>\n",
			no_output_of_prior_pics_flag);
	}
	unsigned int slice_pic_parameter_set_id = get_ue(r->CPB, &shift, 63);
	printf("<li%s>slice_pic_parameter_set_id: <code>%u</code></li>\n",
		red_if(slice_pic_parameter_set_id >= 4), slice_pic_parameter_set_id);
	if (slice_pic_parameter_set_id >= 4 ||
		r->PPSs[slice_pic_parameter_set_id].num_ref_idx_active[0] == 0)
		return;
	Rage265_slice s = {.p = r->PPSs[slice_pic_parameter_set_id]};
	unsigned int dependent_slice_segment_flag = 0;
	if (!first_slice_segment_in_pic_flag) {
		if (s.p.dependent_slice_segments_enabled_flag) {
			dependent_slice_segment_flag = get_u1(r->CPB, &shift);
			printf("<li%s>dependent_slice_segment_flag: <code>%x</code></li>\n",
				red_if(dependent_slice_segment_flag), dependent_slice_segment_flag);
		}
		unsigned int last_ctb = s.p.PicWidthInCtbsY * s.p.PicHeightInCtbsY - 1;
		if (last_ctb > 0) {
			unsigned int slice_segment_address = min(get_uv(r->CPB, &shift,
				WORD_BIT - __builtin_clz(last_ctb)), last_ctb);
			s.ctb_x = slice_segment_address % s.p.PicWidthInCtbsY;
			s.ctb_y = slice_segment_address / s.p.PicWidthInCtbsY;
			printf("<li>slice_segment_address: <code>%u</code></li>\n",
				slice_segment_address);
		}
	}
	if (!dependent_slice_segment_flag) {
		shift += s.p.num_extra_slice_header_bits;
		s.slice_type = get_ue(r->CPB, &shift, 2);
		printf("<li>slice_type: <code>%u (%s)</code></li>\n",
			s.slice_type, slice_type_names[s.slice_type]);
		unsigned int pic_output_flag = 0;
		if (s.p.output_flag_present_flag) {
			pic_output_flag = get_u1(r->CPB, &shift);
			printf("<li>pic_output_flag: <code>%x</code></li>\n", pic_output_flag);
		}
		if (s.p.separate_colour_plane_flag) {
			s.colour_plane_id = get_uv(r->CPB, &shift, 2);
			printf("<li>colour_plane_id: <code>%u (%s)</code></li>\n",
				s.colour_plane_id, colour_plane_id_names[s.colour_plane_id]);
		}
		int PicOrderCntVal = 0;
		unsigned int NumPocTotalCurr = 0;
		if (r->nal_unit_type != 19 && r->nal_unit_type != 20) {
			
			/* 8.3.1 Decoding process for picture order count */
			unsigned int slice_pic_order_cnt_lsb = get_uv(r->CPB, &shift,
				s.p.log2_max_pic_order_cnt_lsb);
			int MaxPicOrderCntLsb = 1 << s.p.log2_max_pic_order_cnt_lsb;
			PicOrderCntVal = (r->prevPicOrderCntVal & -MaxPicOrderCntLsb) |
				slice_pic_order_cnt_lsb;
			if (r->prevPicOrderCntVal - PicOrderCntVal >= MaxPicOrderCntLsb / 2)
				PicOrderCntVal += MaxPicOrderCntLsb;
			if (PicOrderCntVal - r->prevPicOrderCntVal > MaxPicOrderCntLsb / 2)
				PicOrderCntVal -= MaxPicOrderCntLsb;
			printf("<li>PicOrderCntVal: <code>%d</code></li>\n", PicOrderCntVal);
			
			/* 8.3.2 Decoding process for reference picture set */
			unsigned int short_term_ref_pic_set_sps_flag = get_u1(r->CPB, &shift);
			printf("<li>short_term_ref_pic_set_sps_flag: <code>%x</code></li>\n",
				short_term_ref_pic_set_sps_flag);
			unsigned int short_term_ref_pic_set_idx = 0;
			if (!short_term_ref_pic_set_sps_flag) {
				short_term_ref_pic_set_idx = s.p.num_short_term_ref_pic_sets;
				shift = parse_short_term_ref_pic_set(r->short_term_RPSs,
					short_term_ref_pic_set_idx, r->CPB, shift,
					short_term_ref_pic_set_idx, &s.p);
			} else if (s.p.num_short_term_ref_pic_sets > 1) {
				short_term_ref_pic_set_idx = min(get_uv(r->CPB, &shift,
					WORD_BIT - __builtin_clz(s.p.num_short_term_ref_pic_sets - 1)),
					s.p.num_short_term_ref_pic_sets);
				printf("<li>short_term_ref_pic_set_idx: <code>%u</code></li>\n",
					short_term_ref_pic_set_idx);
			}
			const uint16_t *st_RPS = r->short_term_RPSs[short_term_ref_pic_set_idx];
			
			if (s.p.long_term_ref_pics_present_flag) {
				unsigned int num_long_term_sps = 0;
				if (s.p.num_long_term_ref_pics_sps > 0) {
					num_long_term_sps = min(get_ue8(r->CPB, &shift),
						min(s.p.num_long_term_ref_pics_sps,
						s.p.max_dec_pic_buffering - 1 - (st_RPS[15] >> 8)));
					printf("<li>num_long_term_sps: <code>%u</code></li>\n",
						num_long_term_sps);
				}
				unsigned int num_long_term_pics = min(get_ue8(r->CPB, &shift),
					s.p.max_dec_pic_buffering - 1 - (st_RPS[15] >> 8) - num_long_term_sps);
				printf("<li>num_long_term_pics: <code>%u</code></li>\n",
					num_long_term_pics);
				unsigned int DeltaPocMsbCycleLt = 0;
				for (unsigned int i = 0; i < num_long_term_sps + num_long_term_pics; i++) {
					unsigned int pocLt, used_by_curr_pic_lt_flag;
					if (i < num_long_term_sps) {
						unsigned int lt_idx_sps = 0;
						if (s.p.num_long_term_ref_pics_sps > 1) {
							unsigned int lt_idx_sps = min(get_uv(r->CPB, &shift,
								WORD_BIT - __builtin_clz(s.p.num_long_term_ref_pics_sps - 1)),
								s.p.num_long_term_ref_pics_sps - 1);
						}
						pocLt = s.p.lt_ref_pic_poc_lsb_sps[lt_idx_sps];
						used_by_curr_pic_lt_flag = (s.p.used_by_curr_pic_lt_sps_flags >> i) & 1;
					} else {
						if (i == num_long_term_sps)
							DeltaPocMsbCycleLt = 0;
						pocLt = get_uv(r->CPB, &shift, s.p.log2_max_pic_order_cnt_lsb);
						used_by_curr_pic_lt_flag = get_u1(r->CPB, &shift);
					}
					unsigned int delta_poc_msb_present_flag = get_u1(r->CPB, &shift);
					if (delta_poc_msb_present_flag) {
						DeltaPocMsbCycleLt += get_ue64(r->CPB, &shift);
						pocLt += (PicOrderCntVal & -MaxPicOrderCntLsb) -
							(DeltaPocMsbCycleLt << s.p.log2_max_pic_order_cnt_lsb);
					}
					printf("<li>pocLt[%u]: <code>%u, %s</code></li>\n",
						i, pocLt, used_by_curr_pic_lt_flag ? "used" : "follow");
				}
			}
			if (s.p.temporal_mvp_enabled_flag) {
				s.p.temporal_mvp_enabled_flag = get_u1(r->CPB, &shift);
				printf("<li>slice_temporal_mvp_enabled_flag: <code>%x</code></li>\n",
					s.p.temporal_mvp_enabled_flag);
			}
		}
		if (s.p.sample_adaptive_offset_enabled_flag) {
			s.slice_sao_luma_flag = get_u1(r->CPB, &shift);
			s.slice_sao_chroma_flag = get_u1(r->CPB, &shift);
			printf("<li>slice_sao_luma_flag: <code>%x</code></li>\n"
				"<li>slice_sao_chroma_flag: <code>%x</code></li>\n",
				s.slice_sao_luma_flag,
				s.slice_sao_chroma_flag);
		}
		if (s.slice_type <= 1) {
			unsigned int num_ref_idx_active_override_flag = get_u1(r->CPB, &shift);
			for (int l = 0; num_ref_idx_active_override_flag && l < 2 - s.slice_type; l++) {
				s.p.num_ref_idx_active[l] = get_ue(r->CPB, &shift, 14) + 1;
				printf("<li>num_ref_idx_l%u_active: <code>%u</code></li>\n",
					l, s.p.num_ref_idx_active[l]);
			}
			if (s.p.lists_modification_present_flag && NumPocTotalCurr > 1) {
				for (int l = 0; l < 2 - s.slice_type; l++)
					shift = parse_ref_pic_list_modification(r->CPB, shift);
			}
			if (s.slice_type == 0) {
				s.mvd_l1_zero_flag = get_u1(r->CPB, &shift);
				printf("<li>mvd_l1_zero_flag: <code>%x</code></li>\n",
					s.mvd_l1_zero_flag);
			}
			unsigned int cabac_init_flag = 0;
			if (s.p.cabac_init_present_flag) {
				cabac_init_flag = get_u1(r->CPB, &shift);
				printf("<li>cabac_init_flag: <code>%x</code></li>\n",
					cabac_init_flag);
			}
			if (s.p.temporal_mvp_enabled_flag) {
				if (s.slice_type == 0) {
					s.collocated_from_l1_flag = get_u1(r->CPB, &shift) ^ 1;
					printf("<li>collocated_from_l0_flag: <code>%x</code></li>\n",
						s.collocated_from_l1_flag ^ 1);
				}
				if (s.p.num_ref_idx_active[s.collocated_from_l1_flag] > 1) {
					s.collocated_ref_idx = min(get_ue8(r->CPB, &shift),
						s.p.num_ref_idx_active[s.collocated_from_l1_flag] - 1);
					printf("<li>collocated_ref_idx: <code>%u</code></li>\n",
						s.collocated_ref_idx);
				}
			}
			if ((s.p.weighted_pred_flags >> s.slice_type) & 1) {
				unsigned int luma_log2_weight_denom = get_ue(r->CPB, &shift, 7);
				printf("<li>luma_log2_weight_denom: <code>%u</code></li>\n",
					luma_log2_weight_denom);
				if (s.p.ChromaArrayType != 0 || s.p.separate_colour_plane_flag) {
					unsigned int codeNum = get_ue8(r->CPB, &shift);
					unsigned int abs = (codeNum + 1) / 2;
					unsigned int sign = (codeNum % 2) - 1;
					unsigned int chroma_log2_weight_denom = luma_log2_weight_denom +
						(abs ^ sign) - sign;
					if (chroma_log2_weight_denom > 7) // unsigned min
						chroma_log2_weight_denom = 7;
					printf("<li>chroma_log2_weight_denom: <code>%u</code></li>\n",
						chroma_log2_weight_denom);
				}
				for (int l = 0; l < 2 - s.slice_type; l++)
					shift = parse_pred_weight_table(r->CPB, shift);
			}
			s.MaxNumMergeCand = 5 - get_ue(r->CPB, &shift, 4);
			printf("<li>MaxNumMergeCand: <code>%u</code></li>\n",
				s.MaxNumMergeCand);
		}
		unsigned int codeNum = get_ue32(r->CPB, &shift);
		unsigned int abs = (codeNum + 1) / 2;
		unsigned int sign = (codeNum % 2) - 1;
		s.p.Qp = min(max(s.p.Qp + (abs ^ sign) - sign, -s.p.QpBdOffset_Y), 51);
		printf("<li>SliceQp<sub>Y</sub>: <code>%d</code></li>\n", s.p.Qp);
		if (s.p.pps_slice_chroma_qp_offsets_present_flag) {
			s.p.cb_qp_offset = get_se(r->CPB, &shift, -12, 12);
			s.p.cr_qp_offset = get_se(r->CPB, &shift, -12, 12);
			printf("<li>slice_cb_qp_offset: <code>%d</code></li>\n"
				"<li>slice_cr_qp_offset: <code>%d</code></li>\n",
				s.p.cb_qp_offset,
				s.p.cr_qp_offset);
		}
		unsigned int deblocking_filter_override_flag = 0;
		if (s.p.deblocking_filter_override_enabled_flag)
			deblocking_filter_override_flag = get_u1(r->CPB, &shift);
		if (deblocking_filter_override_flag) {
			s.p.deblocking_filter_disabled_flag = get_u1(r->CPB, &shift);
			printf("<li>slice_deblocking_filter_disabled_flag: <code>%x</code></li>\n",
				s.p.deblocking_filter_disabled_flag);
			if (!s.p.deblocking_filter_disabled_flag) {
				s.p.beta_offset = get_se(r->CPB, &shift, -6, 6) * 2;
				s.p.tc_offset = get_se(r->CPB, &shift, -6, 6) * 2;
				printf("<li>slice_beta_offset: <code>%d</code></li>\n"
					"<li>slice_tc_offset: <code>%d</code></li>\n",
					s.p.beta_offset,
					s.p.tc_offset);
			}
		}
		if (s.p.loop_filter_across_slices_enabled_flag && (s.slice_sao_luma_flag ||
			s.slice_sao_chroma_flag || !s.p.deblocking_filter_disabled_flag)) {
			s.p.loop_filter_across_slices_enabled_flag = get_u1(r->CPB, &shift);
			printf("<li>slice_loop_filter_across_slices_enabled_flag: <code>%x</code></li>\n",
				s.p.loop_filter_across_slices_enabled_flag);
		}
	}
	if (s.p.num_tile_columns > 1 || s.p.num_tile_rows > 1 ||
		s.p.entropy_coding_sync_enabled_flag) {
		/* To be refined once multi-threading is implemented. */
		unsigned int num_entry_point_offsets = get_ue(r->CPB, &shift, 23188);
		printf("<li>num_entry_point_offsets: <code>%u</code></li>\n",
			num_entry_point_offsets);
		if (num_entry_point_offsets > 0) {
			unsigned int offset_len = get_ue(r->CPB, &shift, 31) + 1;
			for (unsigned int i = 0; i < num_entry_point_offsets; i++) {
				unsigned int entry_point_offset = get_uv(r->CPB, &shift, offset_len) + 1;
				printf("<li>entry_point_offset[%u]: <code>%u</code></li>\n",
					i, entry_point_offset);
			}
		}
	}
	if (s.p.slice_segment_header_extension_present_flag)
		shift += get_ue(r->CPB, &shift, 256);
	s.c.shift = (shift + 8) & -8;
	if (dependent_slice_segment_flag)
		return;
}



static void parse_AUD(Rage265_ctx *r, unsigned int lim) {
	static const char * const pic_type_names[8] = {"I", "P, I", "B, P, I", [3 ... 7] = "unknown"};
	unsigned int pic_type = *r->CPB >> 5;
	printf("<li%s>pic_type: <code>%u (%s)</code></li>\n",
		red_if(lim != 3), pic_type, pic_type_names[pic_type]);
}



static unsigned int parse_scaling_list_data(Rage265_parameter_set *p, const uint8_t *CPB, unsigned int shift) {
	for (unsigned int matrixId = 0; matrixId < 6; matrixId++) {
		printf("<li>ScalingList[0][%u]: <code>", matrixId);
		unsigned int scaling_list_pred_mode_flag = get_u1(CPB, &shift);
		if (!scaling_list_pred_mode_flag) {
			unsigned int refMatrixId = matrixId - min(get_ue8(CPB, &shift), matrixId);
			if (refMatrixId != matrixId)
				memcpy(p->ScalingFactor4x4[matrixId], p->ScalingFactor4x4[refMatrixId], 16);
			else
				memset(p->ScalingFactor4x4[matrixId], 16, 16);
			const char *str = (refMatrixId != matrixId) ? "ScalingList[0][%u]" : "default";
			printf(str, refMatrixId);
		} else for (unsigned int nextCoef = 8, i = 0; i < 16; i++) {
			nextCoef = (nextCoef + get_se(CPB, &shift, -128, 127)) & 0xff;
			p->ScalingFactor4x4[matrixId][ScanOrder4x4[0][i]] = nextCoef;
			printf(" %u", nextCoef);
		}
		printf("</code></li>\n");
	}
	for (unsigned int sizeId = 0; sizeId < 3; sizeId++) {
		unsigned int num = (sizeId == 2) ? 2 : 6;
		for (unsigned int matrixId = 0; matrixId < num; matrixId++) {
			printf("<li>ScalingList[%u][%u]: <code>", sizeId + 1, matrixId);
			unsigned int scaling_list_pred_mode_flag = get_u1(CPB, &shift);
			if (!scaling_list_pred_mode_flag) {
				unsigned int refMatrixId = matrixId - min(get_ue8(CPB, &shift), matrixId);
				memcpy((&p->ScalingFactor8x8)[sizeId][matrixId],
					(refMatrixId != matrixId) ? (&p->ScalingFactor8x8)[sizeId][refMatrixId] :
					(matrixId < num / 2) ? intra_ScalingFactor : inter_ScalingFactor, 64);
				const char *str = (refMatrixId != matrixId) ? "ScalingList[%u][%u]" : "default";
				printf(str, sizeId + 1, refMatrixId);
			} else {
				unsigned int nextCoef = 8;
				if (sizeId > 0) {
					nextCoef = 8 + get_se(CPB, &shift, -7, 247);
					(sizeId == 1 ? p->ScalingFactor4x4[matrixId] : p->ScalingFactor8x8[matrixId])[0] = nextCoef;
					printf("%u(DC)", nextCoef);
				}
				for (unsigned int i = 0; i < 64; i++) {
					nextCoef = (nextCoef + get_se(CPB, &shift, -128, 127)) & 0xff;
					(&p->ScalingFactor8x8)[sizeId][matrixId][ScanOrder8x8[0][i]] = nextCoef;
					printf(" %u", nextCoef);
				}
			}
			printf("</code></li>\n");
		}
	}
	return shift;
}



/**
 * This function parses the PPS into a copy of the current SPS, and stores it
 * if no error was detected.
 */
static void parse_PPS(Rage265_ctx *r, unsigned int lim) {
	Rage265_parameter_set p = r->SPS;
	unsigned int shift = 0;
	unsigned int pps_pic_parameter_set_id = get_ue(r->CPB, &shift, 63);
	unsigned int pps_seq_parameter_set_id = get_ue(r->CPB, &shift, 15);
	p.dependent_slice_segments_enabled_flag = get_u1(r->CPB, &shift);
	p.output_flag_present_flag = get_u1(r->CPB, &shift);
	p.num_extra_slice_header_bits = get_uv(r->CPB, &shift, 3);
	p.sign_data_hiding_enabled_flag = get_u1(r->CPB, &shift);
	p.cabac_init_present_flag = get_u1(r->CPB, &shift);
	p.num_ref_idx_active[0] = get_ue(r->CPB, &shift, 14) + 1;
	p.num_ref_idx_active[1] = get_ue(r->CPB, &shift, 14) + 1;
	p.Qp = get_se(r->CPB, &shift, -62, 25) + 26;
	p.constrained_intra_pred_flag = get_u1(r->CPB, &shift);
	p.transform_skip_enabled_flag = get_u1(r->CPB, &shift);
	p.cu_qp_delta_enabled_flag = get_u1(r->CPB, &shift);
	printf("<li%s>pps_pic_parameter_set_id: <code>%u</code></li>\n"
		"<li%s>pps_seq_parameter_set_id: <code>%u</code></li>\n"
		"<li>dependent_slice_segments_enabled_flag: <code>%x</code></li>\n"
		"<li>output_flag_present_flag: <code>%x</code></li>\n"
		"<li>num_extra_slice_header_bits: <code>%u</code></li>\n"
		"<li>sign_data_hiding_enabled_flag: <code>%x</code></li>\n"
		"<li>cabac_init_present_flag: <code>%x</code></li>\n"
		"<li>num_ref_idx_l0_default_active: <code>%u</code></li>\n"
		"<li>num_ref_idx_l1_default_active: <code>%u</code></li>\n"
		"<li>init_qp: <code>%d</code></li>\n"
		"<li>constrained_intra_pred_flag: <code>%x</code></li>\n"
		"<li>transform_skip_enabled_flag: <code>%x</code></li>\n"
		"<li>cu_qp_delta_enabled_flag: <code>%x</code></li>\n",
		red_if(pps_pic_parameter_set_id >= 4), pps_pic_parameter_set_id,
		red_if(pps_seq_parameter_set_id > 0), pps_seq_parameter_set_id,
		p.dependent_slice_segments_enabled_flag,
		p.output_flag_present_flag,
		p.num_extra_slice_header_bits,
		p.sign_data_hiding_enabled_flag,
		p.cabac_init_present_flag,
		p.num_ref_idx_active[0],
		p.num_ref_idx_active[1],
		p.Qp,
		p.constrained_intra_pred_flag,
		p.transform_skip_enabled_flag,
		p.cu_qp_delta_enabled_flag);
	if (p.cu_qp_delta_enabled_flag) {
		p.Log2MinCuQpDeltaSize = max(p.CtbLog2SizeY - get_ue8(r->CPB, &shift), p.MinCbLog2SizeY);
		printf("<li>Log2MinCuQpDeltaSize: <code>%u</code></li>\n",
			p.Log2MinCuQpDeltaSize);
	}
	p.cb_qp_offset = get_se(r->CPB, &shift, -12, 12);
	p.cr_qp_offset = get_se(r->CPB, &shift, -12, 12);
	p.pps_slice_chroma_qp_offsets_present_flag = get_u1(r->CPB, &shift);
	p.weighted_pred_flags = get_u1(r->CPB, &shift) << 1;
	p.weighted_pred_flags |= get_u1(r->CPB, &shift);
	p.transquant_bypass_enabled_flag = get_u1(r->CPB, &shift);
	unsigned int tiles_enabled_flag = get_u1(r->CPB, &shift);
	p.entropy_coding_sync_enabled_flag = get_u1(r->CPB, &shift);
	printf("<li>pps_cb_qp_offset: <code>%d</code></li>\n"
		"<li>pps_cr_qp_offset: <code>%d</code></li>\n"
		"<li>pps_slice_chroma_qp_offsets_present_flag: <code>%x</code></li>\n"
		"<li>weighted_pred_flag: <code>%x</code></li>\n"
		"<li>weighted_bipred_flag: <code>%x</code></li>\n"
		"<li>transquant_bypass_enabled_flag: <code>%x</code></li>\n"
		"<li>tiles_enabled_flag: <code>%x</code></li>\n"
		"<li>entropy_coding_sync_enabled_flag: <code>%x</code></li>\n",
		p.cb_qp_offset,
		p.cr_qp_offset,
		p.pps_slice_chroma_qp_offsets_present_flag,
		p.weighted_pred_flags >> 1,
		p.weighted_pred_flags & 1,
		p.transquant_bypass_enabled_flag,
		tiles_enabled_flag,
		p.entropy_coding_sync_enabled_flag);
	p.colBd[1] = p.PicWidthInCtbsY;
	p.rowBd[1] = p.PicHeightInCtbsY;
	if (tiles_enabled_flag) {
		p.num_tile_columns = min(get_ue(r->CPB, &shift, 21) + 1, p.PicWidthInCtbsY);
		p.num_tile_rows = min(get_ue(r->CPB, &shift, 19) + 1, p.PicHeightInCtbsY);
		p.colBd[p.num_tile_columns] = p.PicWidthInCtbsY;
		p.rowBd[p.num_tile_rows] = p.PicHeightInCtbsY;
		unsigned int uniform_spacing_flag = get_u1(r->CPB, &shift);
		printf("<li>num_tile_columns: <code>%u</code></li>\n"
			"<li>num_tile_rows: <code>%u</code></li>\n"
			"<li>uniform_spacing_flag: <code>%x</code></li>\n",
			p.num_tile_columns,
			p.num_tile_rows,
			uniform_spacing_flag);
		if (uniform_spacing_flag) {
			for (unsigned int i = 1; i < p.num_tile_columns; i++)
				p.colBd[i] = i * p.PicWidthInCtbsY / p.num_tile_columns;
			for (unsigned int i = 1; i < p.num_tile_rows; i++)
				p.rowBd[i] = i * p.PicHeightInCtbsY / p.num_tile_rows;
		} else {
			for (unsigned int i = 1; i < p.num_tile_columns; i++) {
				unsigned int column_width = min(get_ue32(r->CPB, &shift) + 1, p.PicWidthInCtbsY - p.colBd[i - 1] - (p.num_tile_columns - i));
				p.colBd[i] = p.colBd[i - 1] + column_width;
				printf("<li>column_width[%u]: <code>%u</code></li>\n", i - 1, column_width);
			}
			for (unsigned int i = 1; i < p.num_tile_rows; i++) {
				unsigned int row_height = min(get_ue32(r->CPB, &shift) + 1, p.PicHeightInCtbsY - p.rowBd[i - 1] - (p.num_tile_rows - i));
				p.rowBd[i] = p.rowBd[i - 1] + row_height;
				printf("<li>row_height[%u]: <code>%u</code></li>\n", i - 1, row_height);
			}
		}
		p.loop_filter_across_tiles_enabled_flag = get_u1(r->CPB, &shift);
		printf("<li>loop_filter_across_tiles_enabled_flag: <code>%x</code></li>\n",
			p.loop_filter_across_tiles_enabled_flag);
	}
	p.loop_filter_across_slices_enabled_flag = get_u1(r->CPB, &shift);
	unsigned int deblocking_filter_control_present_flag = get_u1(r->CPB, &shift);
	printf("<li>pps_loop_filter_across_slices_enabled_flag: <code>%x</code></li>\n",
		p.loop_filter_across_slices_enabled_flag);
	if (deblocking_filter_control_present_flag) {
		p.deblocking_filter_override_enabled_flag = get_u1(r->CPB, &shift);
		p.deblocking_filter_disabled_flag = get_u1(r->CPB, &shift);
		printf("<li>deblocking_filter_override_enabled_flag: <code>%x</code></li>\n"
			"<li>pps_deblocking_filter_disabled_flag: <code>%x</code></li>\n",
			p.deblocking_filter_override_enabled_flag,
			p.deblocking_filter_disabled_flag);
		if (!p.deblocking_filter_disabled_flag) {
			p.beta_offset = get_se(r->CPB, &shift, -6, 6) * 2;
			p.tc_offset = get_se(r->CPB, &shift, -6, 6) * 2;
			printf("<li>pps_beta_offset: <code>%d<code></li>\n"
				"<li>pps_tc_offset: <code>%d</code></li>\n",
				p.beta_offset,
				p.tc_offset);
		}
	}
	unsigned int pps_scaling_list_data_present_flag = get_u1(r->CPB, &shift);
	if (pps_scaling_list_data_present_flag)
		shift = parse_scaling_list_data(&p, r->CPB, shift);
	p.lists_modification_present_flag = get_u1(r->CPB, &shift);
	p.Log2ParMrgLevel = min(get_ue8(r->CPB, &shift) + 2, p.CtbLog2SizeY);
	p.slice_segment_header_extension_present_flag = get_u1(r->CPB, &shift);
	unsigned int pps_extension_flag = get_u1(r->CPB, &shift);
	printf("<li>lists_modification_present_flag: <code>%x</code></li>\n"
		"<li>Log2ParMrgLevel: <code>%u</code></li>\n"
		"<li>slice_segment_header_extension_present_flag: <code>%x</code></li>\n"
		"<li>pps_extension_flag: <code>%x</code></li>\n",
		p.lists_modification_present_flag,
		p.Log2ParMrgLevel,
		p.slice_segment_header_extension_present_flag,
		pps_extension_flag);
	if (pps_extension_flag && shift < lim)
		shift = lim;
	if (shift != lim)
		printf("<li style=\"color: red\">Bitstream overflow (%d bits)</li>\n", shift - lim);
	
	/* The test for seq_parameter_set_id must happen before any use of SPS data. */
	if (shift == lim && pps_pic_parameter_set_id < 4 && pps_seq_parameter_set_id == 0 && r->DPB != NULL)
		r->PPSs[pps_pic_parameter_set_id] = p;
}



static unsigned int parse_vui_parameters(const uint8_t *CPB, unsigned int shift) {
	
	return shift;
}



static const uint8_t *parse_profile_tier_level(Rage265_parameter_set *p, const uint8_t *c) {
	static const char * const general_profile_idc_names[32] = {"unknown", "Main",
		"Main 10", "Main Still Picture", [4 ... 31] = "unknown"};
	
	p->general_profile_space = c[0] >> 6;
	unsigned int general_tier_flag = (c[0] >> 5) & 1;
	unsigned int general_profile_idc = c[0] & 0x1f;
	p->general_progressive_source_flag = c[5] >> 7;
	p->general_interlaced_source_flag = (c[5] >> 6) & 1;
	unsigned int general_non_packed_constraint_flag = (c[5] >> 5) & 1;
	unsigned int general_frame_only_constraint_flag = (c[5] >> 4) & 1;
	unsigned int general_level_idc = c[11];
	printf("<li%s>general_profile_space: <code>%u</code></li>\n"
		"<li>general_tier_flag: <code>%x</code></li>\n"
		"<li>general_profile_idc: <code>%u (%s profile)</code></li>\n"
		"<li>general_progressive_source_flag: <code>%x</code></li>\n"
		"<li>general_interlaced_source_flag: <code>%x</code></li>\n"
		"<li>general_non_packed_constraint_flag: <code>%x</code></li>\n"
		"<li>general_frame_only_constraint_flag: <code>%x</code></li>\n"
		"<li>general_level_idc: <code>%f</code></li>\n",
		red_if(p->general_profile_space > 0), p->general_profile_space,
		general_tier_flag,
		general_profile_idc, general_profile_idc_names[general_profile_idc],
		p->general_progressive_source_flag,
		p->general_interlaced_source_flag,
		general_non_packed_constraint_flag,
		general_frame_only_constraint_flag,
		(float)general_level_idc / 30);
	unsigned int sub_layer_flags = (c[12] << 8) | c[13];
	/* Sub-layer selection must occur at the demux level, hence any such info is ignored. */
	return c + 14 + 11 * __builtin_popcount(sub_layer_flags & 0xaaa0) +
		__builtin_popcount(sub_layer_flags & 0x5550);
}



/**
 * This function parses the SPS into a Rage265_parameter_set structure, and
 * stores it if no error was detected.
 */
static void parse_SPS(Rage265_ctx *r, unsigned int lim) {
	static const char * const chroma_format_idc_names[4] = {"4:0:0", "4:2:0", "4:2:2", "4:4:4"};
	
	Rage265_parameter_set s = {0};
	s.num_tile_columns = 1;
	s.num_tile_rows = 1;
	s.loop_filter_across_tiles_enabled_flag = 1;
	unsigned int sps_video_parameter_set_id = r->CPB[0] >> 4;
	s.max_sub_layers = min((r->CPB[0] >> 1) & 0x7, 6) + 1;
	s.temporal_id_nesting_flag = r->CPB[0] & 1;
	printf("<li>sps_video_parameter_set_id: <code>%u</code></li>\n"
		"<li>sps_max_sub_layers: <code>%u</code></li>\n"
		"<li>sps_temporal_id_nesting_flag: <code>%x</code></li>\n",
		sps_video_parameter_set_id,
		s.max_sub_layers,
		s.temporal_id_nesting_flag);
	unsigned int shift = 8 * (parse_profile_tier_level(&s, r->CPB + 1) - r->CPB);
	unsigned int sps_seq_parameter_set_id = get_ue(r->CPB, &shift, 15);
	s.ChromaArrayType = get_ue(r->CPB, &shift, 3);
	printf("<li%s>sps_seq_parameter_set_id: <code>%u</code></li>\n"
		"<li>chroma_format_idc: <code>%u (%s)</code></li>\n",
		red_if(sps_seq_parameter_set_id > 0), sps_seq_parameter_set_id,
		s.ChromaArrayType, chroma_format_idc_names[s.ChromaArrayType]);
	if (s.ChromaArrayType == 3) {
		s.separate_colour_plane_flag = get_u1(r->CPB, &shift);
		s.ChromaArrayType &= s.separate_colour_plane_flag - 1;
		printf("<li>separate_colour_plane_flag: <code>%x</code></li>\n",
			s.separate_colour_plane_flag);
	}
	s.pic_width_in_luma_samples = get_ue(r->CPB, &shift, 16888);
	s.pic_height_in_luma_samples = get_ue(r->CPB, &shift, 16888);
	unsigned int conformance_window_flag = get_u1(r->CPB, &shift);
	printf("<li>pic_width_in_luma_samples: <code>%u</code></li>\n"
		"<li>pic_height_in_luma_samples: <code>%u</code></li>\n",
		s.pic_width_in_luma_samples,
		s.pic_height_in_luma_samples);
	if (conformance_window_flag) {
		unsigned int shiftX = (s.ChromaArrayType == 1 || s.ChromaArrayType == 2);
		unsigned int shiftY = (s.ChromaArrayType == 1);
		s.conf_win_left_offset = min(get_ue32(r->CPB, &shift) << shiftX, s.pic_width_in_luma_samples);
		s.conf_win_right_offset = min(get_ue32(r->CPB, &shift) << shiftX, s.pic_width_in_luma_samples - s.conf_win_left_offset);
		s.conf_win_top_offset = min(get_ue32(r->CPB, &shift) << shiftY, s.pic_height_in_luma_samples);
		s.conf_win_bottom_offset = min(get_ue32(r->CPB, &shift) << shiftY, s.pic_height_in_luma_samples - s.conf_win_top_offset);
		printf("<li>conf_win_left_offset: <code>%u</code></li>\n"
			"<li>conf_win_right_offset: <code>%u</code></li>\n"
			"<li>conf_win_top_offset: <code>%u</code></li>\n"
			"<li>conf_win_bottom_offset: <code>%u</code></li>\n",
			s.conf_win_left_offset,
			s.conf_win_right_offset,
			s.conf_win_top_offset,
			s.conf_win_bottom_offset);
	}
	s.BitDepth_Y = get_ue(r->CPB, &shift, 6) + 8;
	s.BitDepth_C = get_ue(r->CPB, &shift, 6) + 8;
	s.QpBdOffset_Y = 6 * (s.BitDepth_Y - 8);
	s.QpBdOffset_C = 6 * (s.BitDepth_C - 8);
	s.log2_max_pic_order_cnt_lsb = get_ue(r->CPB, &shift, 12) + 4;
	unsigned int sps_sub_layer_ordering_info_present_flag = get_u1(r->CPB, &shift);
	printf("<li>BitDepth<sub>Y</sub>: <code>%u</code></li>\n"
		"<li>BitDepth<sub>C</sub>: <code>%u</code></li>\n"
		"<li>MaxPicOrderCntLsb: <code>%u</code></li>\n",
		s.BitDepth_Y,
		s.BitDepth_C,
		1 << s.log2_max_pic_order_cnt_lsb);
	for (unsigned int i = (s.max_sub_layers - 1) & -sps_sub_layer_ordering_info_present_flag; i < s.max_sub_layers; i++) {
		s.max_dec_pic_buffering = get_ue(r->CPB, &shift, 15) + 1;
		s.max_num_reorder_pics = get_ue(r->CPB, &shift, 15);
		unsigned int sps_max_latency_increase = get_ue(r->CPB, &shift, 4294967294) - 1;
		printf("<ul>\n"
			"<li>sps_max_dec_pic_buffering[%u]: <code>%u</code></li>\n"
			"<li>sps_max_num_reorder_pics[%u]: <code>%u</code></li>\n"
			"<li>sps_max_latency_increase[%u]: <code>%u</code></li>\n"
			"</ul>\n",
			i, s.max_dec_pic_buffering,
			i, s.max_num_reorder_pics,
			i, sps_max_latency_increase);
	}
	s.MinCbLog2SizeY = get_ue(r->CPB, &shift, 3) + 3;
	s.CtbLog2SizeY = min(max(s.MinCbLog2SizeY + get_ue8(r->CPB, &shift), 4), 6);
	s.PicWidthInCtbsY = s.pic_width_in_luma_samples >> s.CtbLog2SizeY;
	s.PicHeightInCtbsY = s.pic_height_in_luma_samples >> s.CtbLog2SizeY;
	s.Log2MinTrafoSize = min(get_ue8(r->CPB, &shift) + 2, min(s.MinCbLog2SizeY, 5));
	s.Log2MaxTrafoSize = min(s.Log2MinTrafoSize + get_ue8(r->CPB, &shift), min(s.CtbLog2SizeY, 5));
	s.max_transform_hierarchy_depth_inter = min(get_ue8(r->CPB, &shift), s.CtbLog2SizeY - s.Log2MinTrafoSize);
	s.max_transform_hierarchy_depth_intra = min(get_ue8(r->CPB, &shift), s.CtbLog2SizeY - s.Log2MinTrafoSize);
	unsigned int scaling_list_enabled_flag = get_u1(r->CPB, &shift);
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
		unsigned int sps_scaling_list_data_present_flag = get_u1(r->CPB, &shift);
		if (!sps_scaling_list_data_present_flag) {
			memset(s.ScalingFactor4x4, 16, sizeof(s.ScalingFactor4x4));
			for (unsigned int sizeId = 0; sizeId < 3; sizeId++) {
				unsigned int num = (sizeId == 2) ? 2 : 6;
				for (unsigned int matrixId = 0; matrixId < num; matrixId++) {
					memcpy(&(s.ScalingFactor8x8)[sizeId][matrixId], (matrixId < num / 2) ?
						intra_ScalingFactor : inter_ScalingFactor, 64);
				}
			}
		} else {
			shift = parse_scaling_list_data(&s, r->CPB, shift);
		}
	} else {
		memset(s.ScalingFactor4x4, 1, sizeof(s.ScalingFactor4x4));
		memset(s.ScalingFactor8x8, 1, sizeof(s.ScalingFactor8x8));
		memset(s.ScalingFactor16x16, 1, sizeof(s.ScalingFactor16x16));
		memset(s.ScalingFactor32x32, 1, sizeof(s.ScalingFactor32x32));
	}
	s.amp_enabled_flag = get_u1(r->CPB, &shift);
	s.sample_adaptive_offset_enabled_flag = get_u1(r->CPB, &shift);
	unsigned int pcm_enabled_flag = get_u1(r->CPB, &shift);
	printf("<li>amp_enabled_flag: <code>%x</code></li>\n"
		"<li>sample_adaptive_offset_enabled_flag: <code>%x</code></li>\n"
		"<li>pcm_enabled_flag: <code>%x</code></li>\n",
		s.amp_enabled_flag,
		s.sample_adaptive_offset_enabled_flag,
		pcm_enabled_flag);
	if (pcm_enabled_flag) {
		unsigned int u = get_uv(r->CPB, &shift, 8);
		s.PcmBitDepth_Y = min((u >> 4) + 1, s.BitDepth_Y);
		s.PcmBitDepth_C = min((u & 0xf) + 1, s.BitDepth_C);
		s.Log2MinIpcmCbSizeY = min(max(get_ue8(r->CPB, &shift) + 3, s.MinCbLog2SizeY), min(s.CtbLog2SizeY, 5));
		s.Log2MaxIpcmCbSizeY = min(s.Log2MinIpcmCbSizeY + get_ue8(r->CPB, &shift), min(s.CtbLog2SizeY, 5));
		s.pcm_loop_filter_disabled_flag = get_u1(r->CPB, &shift);
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
	s.num_short_term_ref_pic_sets = get_ue(r->CPB, &shift, 64);
	printf("<li>num_short_term_ref_pic_sets: <code>%u</code></li>\n",
		s.num_short_term_ref_pic_sets);
	uint16_t short_term_RPSs[s.num_short_term_ref_pic_sets][16];
	for (unsigned int i = 0; i < s.num_short_term_ref_pic_sets; i++) {
		shift = parse_short_term_ref_pic_set(short_term_RPSs, i, r->CPB, shift,
			s.num_short_term_ref_pic_sets, &s);
	}
	s.long_term_ref_pics_present_flag = get_u1(r->CPB, &shift);
	printf("<li>long_term_ref_pics_present_flag: <code>%x</code></li>\n",
		s.long_term_ref_pics_present_flag);
	if (s.long_term_ref_pics_present_flag) {
		s.num_long_term_ref_pics_sps = get_ue(r->CPB, &shift, 32);
		printf("<li>num_long_term_ref_pics_sps: <code>%u</code></li>\n",
			s.num_long_term_ref_pics_sps);
		for (unsigned int i = 0; i < s.num_long_term_ref_pics_sps; i++) {
			s.lt_ref_pic_poc_lsb_sps[i] = get_uv(r->CPB, &shift, s.log2_max_pic_order_cnt_lsb);
			s.used_by_curr_pic_lt_sps_flags |= get_u1(r->CPB, &shift) << i;
			printf("<ul>\n"
				"<li>lt_ref_pic_poc_lsb_sps[%u]: <code>%u</code></li>\n"
				"<li>used_by_curr_pic_lt_sps_flag[%u]: <code>%u</code></li>\n"
				"</ul>\n",
				i, s.lt_ref_pic_poc_lsb_sps[i],
				i, (s.used_by_curr_pic_lt_sps_flags >> i) & 1);
		}
	}
	s.temporal_mvp_enabled_flag = get_u1(r->CPB, &shift);
	s.strong_intra_smoothing_enabled_flag = get_u1(r->CPB, &shift);
	unsigned int vui_parameters_present_flag = get_u1(r->CPB, &shift);
	printf("<li>sps_temporal_mvp_enabled_flag: <code>%x</code></li>\n"
		"<li>strong_intra_smoothing_enabled_flag: <code>%x</code></li>\n",
		s.temporal_mvp_enabled_flag,
		s.strong_intra_smoothing_enabled_flag);
	if (vui_parameters_present_flag)
		shift = parse_vui_parameters(r->CPB, shift);
	unsigned int sps_extension_flag = get_u1(r->CPB, &shift);
	printf("<li>sps_extension_flag: <code>%x</code></li>\n", sps_extension_flag);
	if (sps_extension_flag && shift < lim)
		shift = lim;
	if (shift != lim)
		printf("<li style=\"color: red\">Bitstream overflow (%d bits)</li>\n", shift - lim);
	
	/* Clear the r->CPB and reallocate the DPB when the image format changes. */
	if (shift != lim || sps_seq_parameter_set_id > 0 || s.general_profile_space > 0)
		return;
	if ((s.ChromaArrayType ^ r->SPS.ChromaArrayType) |
		(s.pic_width_in_luma_samples ^ r->SPS.pic_width_in_luma_samples) |
		(s.pic_height_in_luma_samples ^ r->SPS.pic_height_in_luma_samples) |
		(s.BitDepth_Y ^ r->SPS.BitDepth_Y) | (s.BitDepth_C ^ r->SPS.BitDepth_C) |
		(s.max_dec_pic_buffering ^ r->SPS.max_dec_pic_buffering)) {
		free((uint8_t *)r->CPB);
		r->CPB = NULL;
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
		if (r->DPB != NULL)
			free(r->DPB);
		r->DPB = calloc(s.max_dec_pic_buffering, sizeof(Rage265_picture) +
			luma_size + 2 * chroma_size);
		uint8_t *p = (uint8_t *)(((uintptr_t)(r->DPB + s.max_dec_pic_buffering) + 63) & -64);
		for (unsigned int i = 0; i < s.max_dec_pic_buffering; i++) {
			r->DPB[i].image = p;
			r->DPB[i].PicOrderCntVal = INT32_MIN;
			p += luma_size + 2 * chroma_size;
		}
		memset(r->PPSs, 0, sizeof(r->PPSs));
	}
	r->SPS = s;
	memcpy(r->short_term_RPSs, short_term_RPSs, sizeof(short_term_RPSs));
}



static unsigned int parse_hrd_parameters(const uint8_t *CPB, unsigned int shift, unsigned int cprms_present_flag, unsigned int max_sub_layers) {
	
	return shift;
}



/**
 * This function currently only prints the VPS to stdout.
 */
static void parse_VPS(Rage265_ctx *r, unsigned int lim) {
	unsigned int u = htobe16(*(uint16_t *)r->CPB);
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
	Rage265_parameter_set v;
	unsigned int shift = 8 * (parse_profile_tier_level(&v, r->CPB + 32) - r->CPB);
	unsigned int vps_sub_layer_ordering_info_present_flag = get_u1(r->CPB, &shift);
	for (unsigned int i = (vps_max_sub_layers - 1) & -vps_sub_layer_ordering_info_present_flag; i < vps_max_sub_layers; i++) {
		unsigned int vps_max_dec_pic_buffering = get_ue(r->CPB, &shift, 15) + 1;
		unsigned int vps_max_num_reorder_pics = get_ue(r->CPB, &shift, 15);
		unsigned int vps_max_latency_increase = get_ue(r->CPB, &shift, 4294967294) - 1;
		printf("<ul>\n"
			"<li>vps_max_dec_pic_buffering[%u]: <code>%u</code></li>\n"
			"<li>vps_max_num_reorder_pics[%u]: <code>%u</code></li>\n"
			"<li>vps_max_latency_increase[%u]: <code>%u</code></li>\n"
			"</ul>\n",
			i, vps_max_dec_pic_buffering,
			i, vps_max_num_reorder_pics,
			i, vps_max_latency_increase);
	}
	unsigned int vps_max_layer_id = get_uv(r->CPB, &shift, 6);
	unsigned int vps_num_layer_sets = get_ue(r->CPB, &shift, 1023) + 1;
	printf("<li>vps_max_layer_id: <code>%u</code></li>\n"
		"<li>vps_num_layer_sets: <code>%u</code></li>\n",
		vps_max_layer_id,
		vps_num_layer_sets);
	for (unsigned int i = 1; i < vps_num_layer_sets; i++) {
		for (unsigned int j = 0; j <= vps_max_layer_id; j++) {
			unsigned int layer_id_included_flag = get_u1(r->CPB, &shift);
			printf("<ul><li>layer_id_included_flag[%u][%u]: <code>%x</code></li></ul>\n",
				i, j, layer_id_included_flag);
		}
	}
	unsigned int vps_timing_info_present_flag = get_u1(r->CPB, &shift);
	if (vps_timing_info_present_flag) {
		unsigned int vps_num_units_in_tick = get_uv(r->CPB, &shift, 32);
		unsigned int vps_time_scale = get_uv(r->CPB, &shift, 32);
		unsigned int vps_poc_proportional_to_timing_flag = get_u1(r->CPB, &shift);
		printf("<li>vps_num_units_in_tick: <code>%u</code></li>\n"
			"<li>vps_time_scale: <code>%u</code></li>\n"
			"<li>vps_poc_proportional_to_timing_flag: <code>%x</code></li>\n",
			vps_num_units_in_tick,
			vps_time_scale,
			vps_poc_proportional_to_timing_flag);
		if (vps_poc_proportional_to_timing_flag) {
			unsigned int vps_num_ticks_poc_diff_one = get_ue(r->CPB, &shift, 4294967294) + 1;
			printf("<li>vps_num_ticks_poc_diff_one: <code>%u</code></li>\n",
				vps_num_ticks_poc_diff_one);
		}
		unsigned int vps_num_hrd_parameters = get_ue(r->CPB, &shift, 1024);
		printf("<li>vps_num_hrd_parameters: <code>%u</code></li>\n",
			vps_num_hrd_parameters);
		for (unsigned int i = 0; i < vps_num_hrd_parameters; i++) {
			unsigned int hrd_layer_set_idx = get_ue(r->CPB, &shift, 1023);
			printf("<ul>\n"
				"<li>hrd_layer_set_idx[%u]: <code>%u</code></li>\n",
				i, hrd_layer_set_idx);
			unsigned int cprms_present_flag = 0;
			if (i > 0)
				cprms_present_flag = get_u1(r->CPB, &shift);
			shift = parse_hrd_parameters(r->CPB, shift, cprms_present_flag, vps_max_sub_layers);
			printf("</ul>\n");
		}
	}
	unsigned int vps_extension_flag = get_u1(r->CPB, &shift);
	printf("<li>vps_extension_flag: <code>%x</code></li>\n", vps_extension_flag);
	if (vps_extension_flag && shift < lim)
		shift = lim;
	if (shift != lim)
		printf("<li style=\"color: red\">Bitstream overflow (%d bits)</li>\n", shift - lim);
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
 * Parse a NAL unit, returning the next picture to output or NULL.
 */
const Rage265_picture *Rage265_parse_NAL(Rage265_ctx *r, const uint8_t *buf, size_t len) {
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
	typedef void (*Parser)(Rage265_ctx *, unsigned int);
	static const Parser parse_nal_unit[64] = {
		[32] = parse_VPS,
		[33] = parse_SPS,
		[34] = parse_PPS,
		[35] = parse_AUD,
	};
	
	/* Allocate the CPB. */
	if (len < 2)
		return NULL;
	const unsigned int suffix_size = 128;
	size_t CPB_size = len - 2 + suffix_size;
	if (CPB_size > 800000000 / 8 + suffix_size) { // Level 6.2, High tier
		CPB_size = 800000000 / 8 + suffix_size;
		len = 800000000 / 8;
	}
	if (r->CPB_size < CPB_size) {
		r->CPB_size = CPB_size;
		if (r->CPB != NULL)
			free((uint8_t *)r->CPB);
		r->CPB = malloc(CPB_size);
		if (r->CPB == NULL)
			return NULL;
	}
	
	/* Copy the CPB while removing every emulation_prevention_three_byte. */
	uint8_t *dst = r->CPB;
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
	unsigned int lim = 8 * (dst - r->CPB) + 7 - __builtin_ctz(*dst);
	memset(dst + 1, 0xff, suffix_size);
	
	/* Parse the nal_unit_header() and branch on nal_unit_type. */
	unsigned int nal_unit_header = (buf[0] << 8) | buf[1];
	r->nal_unit_type = nal_unit_header >> 9;
	unsigned int nuh_layer_id = (nal_unit_header >> 3) & 0x3f;
	unsigned int nuh_temporal_id = (nal_unit_header - 1) & 0x7;
	printf("<ul class=\"frame\">\n"
		"<li%s>nal_unit_type: <code>%u (%s)</code></li>\n"
		"<li>nuh_layer_id: <code>%u</code></li>\n"
		"<li>nuh_temporal_id: <code>%u</code></li>\n",
		red_if(parse_nal_unit[r->nal_unit_type] == NULL), r->nal_unit_type, nal_unit_type_names[r->nal_unit_type],
		nuh_layer_id,
		nuh_temporal_id);
	if (nuh_layer_id == 0 && parse_nal_unit[r->nal_unit_type] != NULL)
		parse_nal_unit[r->nal_unit_type](r, lim);
	printf("</ul>\n");
	return NULL;
}
