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
 * Stores one short term reference picture set at position stRpsIdx, and returns
 * the last shift value. Overflows for at most 30 set bits.
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
			RefRpsIdx -= umin(get_raw_ue(CPB, &shift, 63), RefRpsIdx);
		unsigned int NumNegativePics = st_RPS[RefRpsIdx][15] & 0xff;
		unsigned int NumDeltaPocs = st_RPS[RefRpsIdx][15] >> 8;
		unsigned int delta_rps_sign = get_u1(CPB, &shift);
		unsigned int abs_delta_rps = get_ue(CPB, &shift, 32767) + 1;
		int deltaRps = (abs_delta_rps ^ -delta_rps_sign) + delta_rps_sign;
		unsigned int used_by_curr_pic_flags = 0;
		unsigned int use_delta_flags = 0xffff;
		for (unsigned int j = 0; j <= NumDeltaPocs; j++) {
			unsigned int reorder = (j < NumNegativePics) ? NumNegativePics - 1 - j :
				(j < NumDeltaPocs) ? j + 1 : NumNegativePics;
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
			if (((use_delta_flags >> j) & 1) && dPoc != 0) {
				st_RPS[stRpsIdx][i++] = (umin(abs(dPoc) - 1, 32767) << 1) |
					((used_by_curr_pic_flags >> j) & 1);
				if (dPoc < 0)
					k = i;
			}
		}
		st_RPS[stRpsIdx][15] = (i << 8) | k;
	} else {
		unsigned int num_negative_pics = umin(get_raw_ue(CPB, &shift, 15),
			p->max_dec_pic_buffering - 1);
		unsigned int NumDeltaPocs = umin(num_negative_pics + get_raw_ue(CPB, &shift, 15),
			p->max_dec_pic_buffering - 1);
		// shift reaches lim here in the worst case of overflow
		st_RPS[stRpsIdx][15] = (NumDeltaPocs << 8) | num_negative_pics;
		int delta_poc_minus1 = -1;
		for (unsigned int i = 0; i < NumDeltaPocs; i++) {
			if (i == num_negative_pics)
				delta_poc_minus1 = -1;
			delta_poc_minus1 = umin(delta_poc_minus1 + get_raw_ue(CPB, &shift, 32767) + 1, 32767);
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



/**
 * Parses the reference list part of a slice header into s->RefPicList, and
 * returns the last shift value. unavailable receives the highest POC of any
 * unavailable reference picture. used_for_reference is a bitfield set for each
 * DPB entry in Curr or Foll. Overflows for at most 285 set bits.
 */
static unsigned int parse_slice_ref_pic_set(Rage265_slice *s, int *unavailable_poc,
	unsigned int *used_for_reference, unsigned int *NumPocTotalCurr,
	Rage265_ctx *r, unsigned int shift)
{
	unsigned int short_term_ref_pic_set_sps_flag = get_u1(r->CPB, &shift);
	printf("<li>short_term_ref_pic_set_sps_flag: <code>%x</code></li>\n",
		short_term_ref_pic_set_sps_flag);
	unsigned int short_term_ref_pic_set_idx = 0;
	if (!short_term_ref_pic_set_sps_flag) {
		short_term_ref_pic_set_idx = s->p.num_short_term_ref_pic_sets;
		shift = parse_short_term_ref_pic_set(r->short_term_RPSs,
			short_term_ref_pic_set_idx, r->CPB, shift,
			short_term_ref_pic_set_idx, &s->p);
	} else if (s->p.num_short_term_ref_pic_sets > 1) {
		short_term_ref_pic_set_idx = umin(get_uv(r->CPB, &shift,
			WORD_BIT - __builtin_clz(s->p.num_short_term_ref_pic_sets - 1)),
			s->p.num_short_term_ref_pic_sets);
		printf("<li>short_term_ref_pic_set_idx: <code>%u</code></li>\n",
			short_term_ref_pic_set_idx);
	}
	const uint16_t *st_RPS = r->short_term_RPSs[short_term_ref_pic_set_idx];
	unsigned int NumPocStCurrBefore = 0;
	for (unsigned int i = 0; i < (st_RPS[15] >> 8); i++) {
		unsigned int NumNegativePics = st_RPS[15] & 0xff;
		unsigned int idx = (i < NumNegativePics) ? NumNegativePics - 1 - i : i;
		int pocSt = s->currPic.PicOrderCntVal + (i < NumNegativePics ?
			-(st_RPS[idx] >> 1) - 1 : (st_RPS[idx] >> 1) + 1);
		unsigned int j = 0;
		while (j < s->p.max_dec_pic_buffering && r->DPB[j].PicOrderCntVal != pocSt)
			j++;
		if (j == s->p.max_dec_pic_buffering) {
			*unavailable_poc = max(pocSt, *unavailable_poc);
		} else {
			*used_for_reference |= 1 << j;
			if (st_RPS[idx] & 1) {
				s->RefPicList[0][*NumPocTotalCurr++] = &r->DPB[j];
				NumPocStCurrBefore += i < NumNegativePics;
			}
		}
	}
	for (unsigned int i = 0; i < *NumPocTotalCurr; i++) {
		s->RefPicList[1][i] = s->RefPicList[0][(i < NumPocStCurrBefore) ?
			i + NumPocStCurrBefore : i - NumPocStCurrBefore];
	}
		
	if (s->p.long_term_ref_pics_present_flag) {
		unsigned int rem = s->p.max_dec_pic_buffering - 1 - (st_RPS[15] >> 8);
		unsigned int num_long_term_sps = 0;
		if (s->p.num_long_term_ref_pics_sps > 0) {
			num_long_term_sps = umin(get_raw_ue(r->CPB, &shift, 15),
				umin(s->p.num_long_term_ref_pics_sps, rem));
			printf("<li>num_long_term_sps: <code>%u</code></li>\n",
				num_long_term_sps);
		}
		unsigned int num_long_term_pics = umin(get_raw_ue(r->CPB, &shift, 15),
			rem - num_long_term_sps);
		// shift reaches lim here in the worst case of overflow
		printf("<li>num_long_term_pics: <code>%u</code></li>\n",
			num_long_term_pics);
		unsigned int DeltaPocMsbCycleLt = 0;
		for (unsigned int i = 0; i < num_long_term_sps + num_long_term_pics; i++) {
			unsigned int pocLt, used_by_curr_pic_lt_flag;
			if (i < num_long_term_sps) {
				unsigned int lt_idx_sps = 0;
				if (s->p.num_long_term_ref_pics_sps > 1) {
					unsigned int lt_idx_sps = umin(get_uv(r->CPB, &shift,
						WORD_BIT - __builtin_clz(s->p.num_long_term_ref_pics_sps - 1)),
						s->p.num_long_term_ref_pics_sps - 1);
				}
				pocLt = s->p.lt_ref_pic_poc_lsb_sps[lt_idx_sps];
				used_by_curr_pic_lt_flag = (s->p.used_by_curr_pic_lt_sps_flags >> i) & 1;
			} else {
				if (i == num_long_term_sps)
					DeltaPocMsbCycleLt = 0;
				pocLt = get_uv(r->CPB, &shift, s->p.log2_max_pic_order_cnt_lsb);
				used_by_curr_pic_lt_flag = get_u1(r->CPB, &shift);
			}
			if (get_u1(r->CPB, &shift)) {
				DeltaPocMsbCycleLt += get_raw_ue(r->CPB, &shift, 1 << 28);
				pocLt += (s->currPic.PicOrderCntVal & (-1 << s->p.log2_max_pic_order_cnt_lsb)) -
					(DeltaPocMsbCycleLt << s->p.log2_max_pic_order_cnt_lsb);
			}
			printf("<li>pocLt[%u]: <code>%u, %s</code></li>\n",
				i, pocLt, used_by_curr_pic_lt_flag ? "used" : "follow");
			unsigned int j = 0;
			while (j < s->p.max_dec_pic_buffering && r->DPB[j].PicOrderCntVal != pocLt)
				j++;
			if (j == s->p.max_dec_pic_buffering) {
				*unavailable_poc = max(pocLt, *unavailable_poc);
				*used_for_reference |= 1 << j;
				if (used_by_curr_pic_lt_flag) {
					s->RefPicList[0][*NumPocTotalCurr] = &r->DPB[j];
					s->RefPicList[1][*NumPocTotalCurr++] = &r->DPB[j];
				}
			}
		}
	}
	return shift;
}



/**
 * Updates s->RefPicList with the permutations present in the slice header, and
 * returns the last shift value. Overflows for at most 122 set bits.
 */
static unsigned int parse_ref_pic_lists_modification(Rage265_slice *s, const uint8_t *CPB, unsigned int shift) {
	// shift reaches lim here in the worst case of overflow
	for (int l = 0; l < 2 - s->slice_type; l++) {
		if (get_u1(CPB, &shift)) {
			const Rage265_picture *list[s->p.num_ref_idx_active[l]];
			for (unsigned int i = 0; i < s->p.num_ref_idx_active[l]; i++) {
				unsigned int list_entry = umin(get_uv(CPB, &shift,
					WORD_BIT - __builtin_clz(s->p.num_ref_idx_active[l] - 1)),
					s->p.num_ref_idx_active[l]);
				list[i] = s->RefPicList[l][list_entry];
				printf("<li>list_entry_l%u[%u]: <code>%u</code></li>\n",
					l, i, list_entry);
			}
			memcpy(s->RefPicList[l], list, sizeof(list));
		}
	}
	return shift;
}



/**
 * Parses the prediction weights and offsets, and returns the last shift value.
 * Overflows for at most 242 set bits.
 */
static unsigned int parse_pred_weight_table(Rage265_slice *s, const uint8_t *CPB, unsigned int shift) {
	// shift reaches lim here in the worst case of overflow
	s->log2_weight_denom[0] = get_ue(CPB, &shift, 7);
	printf("<li>luma_log2_weight_denom: <code>%u</code></li>\n",
		s->log2_weight_denom[0]);
	if (s->p.chroma_format_idc != 0) {
		s->log2_weight_denom[1] = s->log2_weight_denom[2] =
			umin(s->log2_weight_denom[0] + get_raw_se(CPB, &shift, -7, 7), 7);
		printf("<li>ChromaLog2WeightDenom: <code>%u</code></li>\n",
			s->log2_weight_denom[1]);
	}
	for (int l = 0; l < 2 - s->slice_type; l++) {
		unsigned int luma_weight_flags = get_uv(CPB, &shift, s->p.num_ref_idx_active[l]);
		unsigned int chroma_weight_flags = (s->p.chroma_format_idc != 0) ?
			get_uv(CPB, &shift, s->p.num_ref_idx_active[l]) : 0;
		for (unsigned int i = 0; i < s->p.num_ref_idx_active[l]; i++) {
			if (luma_weight_flags & (1 << (s->p.num_ref_idx_active[l] - 1 - i))) {
				s->delta_weights[l][i][0] = get_se(CPB, &shift, -128, 127);
				s->delta_offsets[l][i][0] = get_se(CPB, &shift, -128, 127);
				printf("<li>delta_luma_weight_l%u[%u]: <code>%d</code></li>\n"
					"<li>luma_offset_l%u[%u]: <code>%d</code></li>\n",
					l, i, s->delta_weights[l][i][0],
					l, i, s->delta_offsets[l][i][0]);
			}
			if (chroma_weight_flags & (1 << (s->p.num_ref_idx_active[l] - 1 - i))) {
				for (unsigned int j = 0; j < 2; j++) {
					s->delta_weights[l][i][1 + j] = get_se(CPB, &shift, -128, 127);
					s->delta_offsets[l][i][1 + j] = min(max(get_raw_se(CPB, &shift, -512, 511) -
						(s->delta_weights[l][i][1 + j] << (7 - s->log2_weight_denom[1])), -128), 127);
					printf("<li>delta_chroma_weight_l%u[%u][%u]: <code>%d</code></li>\n"
						"<li>ChromaOffsetL%u[%u][%u]: <code>%d</code></li>\n",
						l, i, j, s->delta_weights[l][i][1 + j],
						l, i, j, s->delta_offsets[l][i][1 + j]);
				}
			}
		}
	}
	return shift;
}



/** Fills an array of samples with 1<<(BitDepth-1). */
static void paint_grey(uint16_t *buf, unsigned int num, unsigned int BitDepth) {
	// 32 is the size of a 4:2:0 chroma block when MinCbLog2SizeY==3
	typedef union { uint16_t h[16]; } Vect __attribute__((aligned(32)));
	Vect v = (Vect){.h = {[0 ... 15] = 1 << (BitDepth - 1)}};
	for (unsigned int u = 0; u < num; u += 16)
		*(Vect *)(buf + u) = v;
}



/**
 * Copies a parameter set and parses a slice header into a Rage265_slice
 * context, and yields decoding to parse_slice_segment_data(). The next picture
 * to be output is returned. Overflows for at most 695 set bits.
 */
static const Rage265_picture *parse_slice_segment_layer(Rage265_ctx *r, unsigned int lim) {
	static const char * const slice_type_names[3] = {"B", "P", "I"};
	static const char * const colour_plane_id_names[3] = {"Y", "Cb", "Cr"};
	
	unsigned int shift = 0;
	// shift reaches lim here in the worst case of overflow
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
		return NULL;
	Rage265_slice s = {.p = r->PPSs[slice_pic_parameter_set_id],
		.currPic.needed_for_output = 1};
	unsigned int dependent_slice_segment_flag = 0;
	if (!first_slice_segment_in_pic_flag) {
		if (s.p.dependent_slice_segments_enabled_flag) {
			dependent_slice_segment_flag = get_u1(r->CPB, &shift);
			printf("<li%s>dependent_slice_segment_flag: <code>%x</code></li>\n",
				red_if(dependent_slice_segment_flag), dependent_slice_segment_flag);
		}
		unsigned int last_ctb = s.p.PicWidthInCtbsY * s.p.PicHeightInCtbsY - 1;
		if (last_ctb > 0) {
			unsigned int slice_segment_address = umin(get_uv(r->CPB, &shift,
				WORD_BIT - __builtin_clz(last_ctb)), last_ctb);
			s.ctb_x = slice_segment_address % s.p.PicWidthInCtbsY;
			s.ctb_y = slice_segment_address / s.p.PicWidthInCtbsY;
			printf("<li>slice_segment_address: <code>%u</code></li>\n",
				slice_segment_address);
		}
	}
	int unavailable_poc = INT32_MIN;
	unsigned int used_for_reference = 0;
	if (!dependent_slice_segment_flag) {
		shift += s.p.num_extra_slice_header_bits;
		s.slice_type = get_ue(r->CPB, &shift, 2);
		printf("<li>slice_type: <code>%u (%s)</code></li>\n",
			s.slice_type, slice_type_names[s.slice_type]);
		if (s.p.output_flag_present_flag) {
			s.currPic.needed_for_output = get_u1(r->CPB, &shift);
			printf("<li>pic_output_flag: <code>%x</code></li>\n",
				s.currPic.needed_for_output);
		}
		if (s.p.separate_colour_plane_flag) {
			s.colour_plane_id = get_uv(r->CPB, &shift, 2);
			printf("<li>colour_plane_id: <code>%u (%s)</code></li>\n",
				s.colour_plane_id, colour_plane_id_names[s.colour_plane_id]);
		}
		unsigned int NumPocTotalCurr = 0;
		if (r->nal_unit_type != 19 && r->nal_unit_type != 20) {
			
			/* 8.3.1 Decoding process for picture order count */
			unsigned int slice_pic_order_cnt_lsb = get_uv(r->CPB, &shift,
				s.p.log2_max_pic_order_cnt_lsb);
			int MaxPicOrderCntLsb = 1 << s.p.log2_max_pic_order_cnt_lsb;
			s.currPic.PicOrderCntVal = (r->prevPicOrderCntVal & -MaxPicOrderCntLsb) |
				slice_pic_order_cnt_lsb;
			if (r->prevPicOrderCntVal - s.currPic.PicOrderCntVal >= MaxPicOrderCntLsb / 2)
				s.currPic.PicOrderCntVal += MaxPicOrderCntLsb;
			if (s.currPic.PicOrderCntVal - r->prevPicOrderCntVal > MaxPicOrderCntLsb / 2)
				s.currPic.PicOrderCntVal -= MaxPicOrderCntLsb;
			printf("<li>PicOrderCntVal: <code>%d</code></li>\n", s.currPic.PicOrderCntVal);
			
			shift = parse_slice_ref_pic_set(&s, &unavailable_poc, &used_for_reference,
				&NumPocTotalCurr, r, shift);
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
			if (s.p.lists_modification_present_flag && NumPocTotalCurr > 1)
				shift = parse_ref_pic_lists_modification(&s, r->CPB, shift);
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
					s.collocated_ref_idx = umin(get_raw_ue(r->CPB, &shift, 14),
						s.p.num_ref_idx_active[s.collocated_from_l1_flag] - 1);
					printf("<li>collocated_ref_idx: <code>%u</code></li>\n",
						s.collocated_ref_idx);
				}
			}
			if ((s.p.weighted_pred_flags >> s.slice_type) & 1)
				shift = parse_pred_weight_table(&s, r->CPB, shift);
			s.MaxNumMergeCand = 5 - get_ue(r->CPB, &shift, 4);
			printf("<li>MaxNumMergeCand: <code>%u</code></li>\n",
				s.MaxNumMergeCand);
		}
		int slice_qp_delta = get_raw_se(r->CPB, &shift, -87, 87);
		s.p.Qp = min(max(s.p.Qp + slice_qp_delta, -s.p.QpBdOffset_Y), 51);
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
		int num_entry_point_offsets = get_ue(r->CPB, &shift, 23188);
		printf("<li>num_entry_point_offsets: <code>%u</code></li>\n",
			num_entry_point_offsets);
		if (num_entry_point_offsets > 0) {
			unsigned int offset_len = get_ue(r->CPB, &shift, 31) + 1;
			// This would require a huge safety suffix without this protection
			num_entry_point_offsets = min(num_entry_point_offsets,
				(int)(lim - shift) >> offset_len);
			for (unsigned int i = 0; i < num_entry_point_offsets; i++) {
				unsigned int entry_point_offset = get_uv(r->CPB, &shift, offset_len) + 1;
				printf("<li>entry_point_offset[%u]: <code>%u</code></li>\n",
					i, entry_point_offset);
			}
		}
	}
	if (s.p.slice_segment_header_extension_present_flag)
		shift += get_ue(r->CPB, &shift, 256);
	if (((r->CPB[shift / 8] << (shift % 8)) & 0xff) != 0x80)
		printf("<li style=\"color: red\">Slice header overflow</li>\n");
	
	/* Proceed to DPB update when we are assured the slice header was correct. */
	if (shift >= lim || ((r->CPB[shift / 8] << (shift % 8)) & 0xff) != 0x80 ||
		dependent_slice_segment_flag)
		return NULL;
	s.c.shift = (shift + 8) & -8;
	if ((1 << r->nal_unit_type) & 0xffa82a && r->TemporalId == 0)
		r->prevPicOrderCntVal = s.currPic.PicOrderCntVal;
	
	/* When some references are unavailable, create a grey picture for replacement. */
	if (unavailable_poc > INT32_MIN) {
		unsigned int i, avail = 0;
		for (i = 0; i < s.p.max_dec_pic_buffering && !r->DPB[i].grey_picture; i++) {
			if (!r->DPB[i].needed_for_output && !(used_for_reference & (1 << i)))
				avail = i;
		}
		if (i == s.p.max_dec_pic_buffering) {
			i = avail;
			r->DPB[i].needed_for_output = 0;
			r->DPB[i].grey_picture = 1;
			r->DPB[i].PicOrderCntVal = unavailable_poc;
			size_t ctb_size = s.p.PicWidthInCtbsY * s.p.PicHeightInCtbsY * sizeof(Rage265_ctb);
			memset(r->DPB[i].CTBs, 0, ctb_size);
			paint_grey(r->DPB[i].image, s.p.image_offsets[1], s.p.BitDepth_Y);
			paint_grey(r->DPB[i].image + s.p.image_offsets[1], s.p.image_offsets[3] -
				s.p.image_offsets[1], s.p.BitDepth_C);
		}
		used_for_reference |= 1 << i;
		for (int l = 0; l < 2 - s.slice_type; l++) {
			for (unsigned int j = 0; j < s.p.num_ref_idx_active[l]; j++) {
				if (s.RefPicList[l][j] == NULL)
					s.RefPicList[l][j] = &r->DPB[i];
			}
		}
	}
	
	/* Select an DPB emplacement for the new picture. */
	Rage265_picture *avail = r->DPB;
	for (unsigned int i = 0; i < s.p.max_dec_pic_buffering; i++) {
		if (!r->DPB[i].needed_for_output && !(used_for_reference & (1 << i)))
			avail = &r->DPB[i];
	}
	s.currPic.image = avail->image;
	s.currPic.CTBs = avail->CTBs;
	*avail = s.currPic;
	
	// TODO: Parse the slice data
	
	/* Select an image for output. */
	Rage265_picture *output = NULL;
	unsigned int num_reordered_pics = 0;
	for (unsigned int i = 0; i < s.p.max_dec_pic_buffering; i++) {
		if (r->DPB[i].needed_for_output) {
			num_reordered_pics++;
			if (output == NULL || r->DPB[i].PicOrderCntVal < output->PicOrderCntVal)
				output = &r->DPB[i];
		}
	}
	if (num_reordered_pics <= s.p.max_num_reorder_pics && r->nal_unit_type < 20 &&
		r->nal_unit_type != 18)
		output = NULL;
	else
		output->needed_for_output = 0;
	return output;
}



/** Frees all internal structures and clears the Rage265_ctx. */
static const Rage265_picture *parse_end_of_bitstream(Rage265_ctx *r, unsigned int lim) {
	if (lim == 0) {
		if (r->CPB != NULL)
			free(r->CPB);
		if (r->DPB[0].image != NULL)
			free(r->DPB[0].image);
		memset(r, 0, sizeof(*r));
	}
	return NULL;
}



/** Returns the next output picture if there is one, otherwise resets prevPicOrderCntVal. */
static const Rage265_picture *parse_end_of_seq(Rage265_ctx *r, unsigned int lim) {
	Rage265_picture *output = NULL;
	if (lim == 0) {
		for (unsigned int i = 0; i < r->SPS.max_dec_pic_buffering; i++) {
			if (r->DPB[i].needed_for_output && (output == NULL ||
				r->DPB[i].PicOrderCntVal < output->PicOrderCntVal))
				output = &r->DPB[i];
		}
		if (output != NULL)
			output->needed_for_output = 0;
		else
			r->prevPicOrderCntVal = 0;
	}
	return output;
}



/** Prints out an AUD. */
static const Rage265_picture *parse_access_unit_delimiter(Rage265_ctx *r, unsigned int lim) {
	static const char * const pic_type_names[8] = {"I", "P, I", "B, P, I", [3 ... 7] = "unknown"};
	unsigned int pic_type = *r->CPB >> 5;
	printf("<li%s>pic_type: <code>%u (%s)</code></li>\n",
		red_if(lim != 3), pic_type, pic_type_names[pic_type]);
	return NULL;
}



/**
 * Parses the scaling lists into p->ScalingFactorNxN, and returns the last shift
 * value. Overflows for at most 1000 set bits.
 */
static unsigned int parse_scaling_list_data(Rage265_parameter_set *p,
	const uint8_t *CPB, unsigned int shift)
{
	// shift reaches lim here in the worst case of overflow
	for (unsigned int matrixId = 0; matrixId < 6; matrixId++) {
		printf("<li>ScalingList[0][%u]: <code>", matrixId);
		if (!get_u1(CPB, &shift)) {
			unsigned int refMatrixId = matrixId - umin(get_raw_ue(CPB, &shift, 6), matrixId);
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
			if (!get_u1(CPB, &shift)) {
				unsigned int refMatrixId = matrixId - umin(get_raw_ue(CPB, &shift, 6), matrixId);
				// The address trick assumes no variable is inserted between the arrays
				memcpy((&p->ScalingFactor8x8)[sizeId][matrixId],
					(refMatrixId != matrixId) ? (&p->ScalingFactor8x8)[sizeId][refMatrixId] :
					(matrixId < num / 2) ? intra_ScalingFactor : inter_ScalingFactor, 64);
				const char *str = (refMatrixId != matrixId) ? "ScalingList[%u][%u]" : "default";
				printf(str, sizeId + 1, refMatrixId);
			} else {
				unsigned int nextCoef = 8;
				if (sizeId > 0) {
					nextCoef = 8 + get_se(CPB, &shift, -7, 247);
					(sizeId == 1 ? p->ScalingFactor4x4[matrixId] :
						p->ScalingFactor8x8[matrixId])[0] = nextCoef;
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
 * Parses the PPS into a copy of the current SPS, and returns NULL.
 * Overflows for at most 1050 set bits.
 */
static const Rage265_picture *parse_picture_parameter_set(Rage265_ctx *r, unsigned int lim) {
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
		p.Log2MinCuQpDeltaSize = umax(p.CtbLog2SizeY - get_raw_ue(r->CPB, &shift, 3),
			p.MinCbLog2SizeY);
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
		p.num_tile_columns = umin(get_ue(r->CPB, &shift, 21) + 1, p.PicWidthInCtbsY);
		p.num_tile_rows = umin(get_ue(r->CPB, &shift, 19) + 1, p.PicHeightInCtbsY);
		p.colBd[p.num_tile_columns] = p.PicWidthInCtbsY;
		p.rowBd[p.num_tile_rows] = p.PicHeightInCtbsY;
		unsigned int uniform_spacing_flag = get_u1(r->CPB, &shift);
		// shift reaches lim here in the worst case of overflow
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
				unsigned int column_width = umin(get_raw_ue(r->CPB, &shift, 1055) + 1,
					p.PicWidthInCtbsY - p.colBd[i - 1] - (p.num_tile_columns - i));
				p.colBd[i] = p.colBd[i - 1] + column_width;
				printf("<li>column_width[%u]: <code>%u</code></li>\n", i - 1, column_width);
			}
			for (unsigned int i = 1; i < p.num_tile_rows; i++) {
				unsigned int row_height = umin(get_raw_ue(r->CPB, &shift, 1055) + 1,
					p.PicHeightInCtbsY - p.rowBd[i - 1] - (p.num_tile_rows - i));
				p.rowBd[i] = p.rowBd[i - 1] + row_height;
				printf("<li>row_height[%u]: <code>%u</code></li>\n", i - 1, row_height);
			}
		}
		p.loop_filter_across_tiles_enabled_flag = get_u1(r->CPB, &shift);
		printf("<li>loop_filter_across_tiles_enabled_flag: <code>%x</code></li>\n",
			p.loop_filter_across_tiles_enabled_flag);
	}
	p.loop_filter_across_slices_enabled_flag = get_u1(r->CPB, &shift);
	printf("<li>pps_loop_filter_across_slices_enabled_flag: <code>%x</code></li>\n",
		p.loop_filter_across_slices_enabled_flag);
	if (get_u1(r->CPB, &shift)) {
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
	if (get_u1(r->CPB, &shift))
		shift = parse_scaling_list_data(&p, r->CPB, shift);
	p.lists_modification_present_flag = get_u1(r->CPB, &shift);
	p.Log2ParMrgLevel = umin(get_raw_ue(r->CPB, &shift, 4) + 2, p.CtbLog2SizeY);
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
	if (shift == lim && pps_pic_parameter_set_id < 4 && pps_seq_parameter_set_id == 0 &&
		r->DPB[0].image != NULL)
		r->PPSs[pps_pic_parameter_set_id] = p;
	return NULL;
}



/** Returns the last shift value. Overflows for at most CpbCnt*5 set bits. */
static unsigned int parse_sub_layer_hrd_parameters(const uint8_t *CPB,
	unsigned int shift, unsigned int CpbCnt, unsigned int sub_pic_hrd_params_present_flag)
{
	// shift reaches lim here in the worst case of overflow
	for (unsigned int i = 0; i < CpbCnt; i++) {
		unsigned int bit_rate_value = get_ue(CPB, &shift, 4294967294) + 1;
		unsigned int cpb_size_value = get_ue(CPB, &shift, 4294967294) + 1;
		printf("<ul>\n"
			"<li>bit_rate_value[%u]: <code>%u</code></li>\n"
			"<li>cpb_size_value[%u]: <code>%u</code></li>\n",
			i, bit_rate_value,
			i, cpb_size_value);
		if (sub_pic_hrd_params_present_flag) {
			unsigned int cpb_size_du_value = get_ue(CPB, &shift, 4294967294) + 1;
			unsigned int bit_rate_du_value = get_ue(CPB, &shift, 4294967294) + 1;
			printf("<li>cpb_size_du_value[%u]: <code>%u</code></li>\n"
				"<li>bit_rate_du_value[%u]: <code>%u</code></li>\n",
				i, cpb_size_du_value,
				i, bit_rate_du_value);
		}
		unsigned int cbr_flag = get_u1(CPB, &shift);
		printf("<li>cbr_flag[%u]: <code>%x</code></li>\n"
			"</ul>\n",
			i, cbr_flag);
	}
	return shift;
}



/** Returns the last shift value. Overflows for at most 398 set bits. */
static unsigned int parse_hrd_parameters(const uint8_t *CPB, unsigned int shift,
	unsigned int cprms_present_flag, unsigned int max_sub_layers)
{
	unsigned int nal_hrd_parameters_present_flag = 0;
	unsigned int vcl_hrd_parameters_present_flag = 0;
	unsigned int sub_pic_hrd_params_present_flag = 0;
	if (cprms_present_flag) {
		nal_hrd_parameters_present_flag = get_u1(CPB, &shift);
		vcl_hrd_parameters_present_flag = get_u1(CPB, &shift);
		if (nal_hrd_parameters_present_flag || vcl_hrd_parameters_present_flag) {
			sub_pic_hrd_params_present_flag = get_u1(CPB, &shift);
			if (sub_pic_hrd_params_present_flag) {
				unsigned int u = get_uv(CPB, &shift, 19);
				unsigned int tick_divisor = (u >> 11) + 2;
				unsigned int du_cpb_removal_delay_increment_length = ((u >> 6) & 0x1f) + 1;
				unsigned int sub_pic_cpb_params_in_pic_timing_sei_flag = (u >> 5) & 1;
				unsigned int dpb_output_delay_du_length = (u & 0x1f) + 1;
				printf("<li>tick_divisor: <code>%u</code></li>\n"
					"<li>du_cpb_removal_delay_increment_length: <code>%u</code></li>\n"
					"<li>sub_pic_cpb_params_in_pic_timing_sei_flag: <code>%x</code></li>\n"
					"<li>dpb_output_delay_du_length: <code>%u</code></li>\n",
					tick_divisor,
					du_cpb_removal_delay_increment_length,
					sub_pic_cpb_params_in_pic_timing_sei_flag,
					dpb_output_delay_du_length);
			}
			unsigned int scale = get_uv(CPB, &shift, 8);
			unsigned int bit_rate_scale = scale >> 4;
			unsigned int cpb_size_scale = scale & 0xf;
			printf("<li>bit_rate_scale: <code>%u</code></li>\n"
				"<li>cpb_size_scale: <code>%u</code></li>\n",
				bit_rate_scale,
				cpb_size_scale);
			if (sub_pic_hrd_params_present_flag) {
				unsigned int cpb_size_du_scale = get_uv(CPB, &shift, 4);
				printf("<li>cpb_size_du_scale: <code>%u</code></li>\n",
					cpb_size_du_scale);
			}
			unsigned int u = get_uv(CPB, &shift, 15);
			unsigned int initial_cpb_removal_delay_length = (u >> 10) + 1;
			unsigned int au_cpb_removal_delay_length = ((u >> 5) & 0x1f) + 1;
			unsigned int dpb_output_delay_length = (u & 0x1f) + 1;
			printf("<li>initial_cpb_removal_delay_length: <code>%u</code></li>\n"
				"<li>au_cpb_removal_delay_length: <code>%u</code></li>\n"
				"dpb_output_delay_length: <code>%u</code></li>\n",
				initial_cpb_removal_delay_length,
				au_cpb_removal_delay_length,
				dpb_output_delay_length);
		}
	}
	for (unsigned int i = 0; i <= max_sub_layers; i++) {
		unsigned int fixed_pic_rate_general_flag = get_u1(CPB, &shift);
		printf("<ul>\n"
			"<li>fixed_pic_rate_general_flag[%u]: <code>%x</code></li>\n",
			i, fixed_pic_rate_general_flag);
		unsigned int fixed_pic_rate_within_cvs_flag = 1;
		if (!fixed_pic_rate_general_flag) {
			fixed_pic_rate_within_cvs_flag = get_u1(CPB, &shift);
			printf("<li>fixed_pic_rate_within_cvs_flag[%u]: <code>%x</code></li>\n",
				i, fixed_pic_rate_within_cvs_flag);
		}
		unsigned int low_delay_hrd_flag = 0;
		if (fixed_pic_rate_within_cvs_flag) {
			unsigned int elemental_duration_in_tc = get_ue(CPB, &shift, 2047) + 1;
			printf("<li>elemental_duration_in_tc[%u]: <code>%u</code></li>\n",
				i, elemental_duration_in_tc);
		} else {
			low_delay_hrd_flag = get_u1(CPB, &shift);
			printf("<li>low_delay_hrd_flag[%u]: <code>%x</code></li>\n",
				i, low_delay_hrd_flag);
		}
		unsigned int CpbCnt = 1;
		if (!low_delay_hrd_flag) {
			CpbCnt = get_ue(CPB, &shift, 31) + 1;
			printf("<li>cpb_cnt[%u]: <code>%u</code></li>\n", i, CpbCnt);
		}
		// shift reaches lim here in the worst case of overflow
		if (nal_hrd_parameters_present_flag)
			shift = parse_sub_layer_hrd_parameters(CPB, shift, CpbCnt, sub_pic_hrd_params_present_flag);
		if (vcl_hrd_parameters_present_flag)
			shift = parse_sub_layer_hrd_parameters(CPB, shift, CpbCnt, sub_pic_hrd_params_present_flag);
		printf("</ul>\n");
	}
	return shift;
}



/** Returns the last shift value. Overflows for at most 559 set bits. */
static unsigned int parse_vui_parameters(const uint8_t *CPB, unsigned int shift,
	unsigned int max_sub_layers)
{
	static const unsigned int ratio2sar[256] = {0, 0x00010001, 0x000c000b,
		0x000a000b, 0x0010000b, 0x00280021, 0x0018000b, 0x0014000b, 0x0020000b,
		0x00500021, 0x0012000b, 0x000f000b, 0x00400021, 0x00a00063, 0x00040003,
		0x00030002, 0x00020001};
	static const char * const video_format_names[8] = {"Component", "PAL",
		"NTSC", "SECAM", "MAC", [5 ... 7] = "Unspecified"};
	
	// shift reaches lim here in the worst case of overflow
	if (get_u1(CPB, &shift)) {
		unsigned int aspect_ratio_idc = get_uv(CPB, &shift, 8);
		unsigned int sar = ratio2sar[aspect_ratio_idc];
		if (aspect_ratio_idc == 255)
			sar = get_uv(CPB, &shift, 32);
		unsigned int sar_width = sar >> 16;
		unsigned int sar_height = sar & 0xffff;
		printf("<li>aspect_ratio: <code>%u (%u:%u)</code></li>\n",
			aspect_ratio_idc, sar_width, sar_height);
	}
	if (get_u1(CPB, &shift)) {
		unsigned int overscan_appropriate_flag = get_u1(CPB, &shift);
		printf("<li>overscan_appropriate_flag: <code>%x</code></li>\n",
			overscan_appropriate_flag);
	}
	if (get_u1(CPB, &shift)) {
		unsigned int video_format = get_uv(CPB, &shift, 3);
		unsigned int video_full_range_flag = get_u1(CPB, &shift);
		printf("<li>video_format: <code>%u (%s)</code></li>\n"
			"<li>video_full_range_flag: <code>%x</code></li>\n",
			video_format, video_format_names[video_format],
			video_full_range_flag);
		if (get_u1(CPB, &shift)) {
			unsigned int u = get_uv(CPB, &shift, 24);
			unsigned int colour_primaries = u >> 16;
			unsigned int transfer_characteristics = (u >> 8) & 0xff;
			unsigned int matrix_coeffs = u & 0xff;
			printf("<li>colour_primaries: <code>%u</code></li>\n"
				"<li>transfer_characteristics: <code>%u</code></li>\n"
				"<li>matrix_coeffs: <code>%u</code></li>\n",
				colour_primaries,
				transfer_characteristics,
				matrix_coeffs);
		}
	}
	if (get_u1(CPB, &shift)) {
		unsigned int chroma_sample_loc_type_top_field = get_ue(CPB, &shift, 5);
		unsigned int chroma_sample_loc_type_bottom_field = get_ue(CPB, &shift, 5);
		printf("<li>chroma_sample_loc_type_top_field: <code>%u</code></li>\n"
			"<li>chroma_sample_loc_type_bottom_field: <code>%u</code></li>\n",
			chroma_sample_loc_type_top_field,
			chroma_sample_loc_type_bottom_field);
	}
	unsigned int neutral_chroma_indication_flag = get_u1(CPB, &shift);
	unsigned int field_seq_flag = get_u1(CPB, &shift);
	unsigned int frame_field_info_present_flag = get_u1(CPB, &shift);
	unsigned int default_display_window_flag = get_u1(CPB, &shift);
	printf("<li>neutral_chroma_indication_flag: <code>%x</code></li>\n"
		"<li>field_seq_flag: <code>%x</code></li>\n"
		"<li>frame_field_info_present_flag: <code>%x</code></li>\n"
		"<li>default_display_window_flag: <code>%x</code></li>\n",
		neutral_chroma_indication_flag,
		field_seq_flag,
		frame_field_info_present_flag,
		default_display_window_flag);
	if (default_display_window_flag) {
		unsigned int def_disp_win_left_offset = get_ue(CPB, &shift, 16887);
		unsigned int def_disp_win_right_offset = get_ue(CPB, &shift, 16887);
		unsigned int def_disp_win_top_offset = get_ue(CPB, &shift, 16887);
		unsigned int def_disp_win_bottom_offset = get_ue(CPB, &shift, 16887);
		printf("<li>def_disp_win_left_offset: <code>%u</code></li>\n"
			"<li>def_disp_win_right_offset: <code>%u</code></li>\n"
			"<li>def_disp_win_top_offset: <code>%u</code></li>\n"
			"<li>def_disp_win_bottom_offset: <code>%u</code></li>\n",
			def_disp_win_left_offset,
			def_disp_win_right_offset,
			def_disp_win_top_offset,
			def_disp_win_bottom_offset);
	}
	if (get_u1(CPB, &shift)) {
		unsigned int vui_num_units_in_tick = get_uv(CPB, &shift, 32);
		unsigned int vui_time_scale = get_uv(CPB, &shift, 32);
		unsigned int vui_poc_proportional_to_timing_flag = get_u1(CPB, &shift);
		printf("<li>vui_num_units_in_tick: <code>%u</code></li>\n"
			"<li>vui_time_scale: <code>%u</code></li>\n"
			"<li>vui_poc_proportional_to_timing_flag: <code>%x</code></li>\n",
			vui_num_units_in_tick,
			vui_time_scale,
			vui_poc_proportional_to_timing_flag);
		if (vui_poc_proportional_to_timing_flag) {
			unsigned int vui_num_ticks_poc_diff_one = get_ue(CPB, &shift, 4294967294) + 1;
			printf("<li>vui_num_ticks_poc_diff_one: <code>%u</code></li>\n",
				vui_num_ticks_poc_diff_one);
		}
		if (get_u1(CPB, &shift))
			shift = parse_hrd_parameters(CPB, shift, 1, max_sub_layers);
	}
	if (get_u1(CPB, &shift)){
		unsigned int tiles_fixed_structure_flag = get_u1(CPB, &shift);
		unsigned int motion_vectors_over_pic_boundaries_flag = get_u1(CPB, &shift);
		unsigned int restricted_ref_pic_lists_flag = get_u1(CPB, &shift);
		unsigned int min_spatial_segmentation_idc = get_ue(CPB, &shift, 4095);
		unsigned int max_bytes_per_pic_denom = get_ue(CPB, &shift, 16);
		unsigned int max_bits_per_min_cu_denom = get_ue(CPB, &shift, 16);
		unsigned int log2_max_mv_length_horizontal = get_ue(CPB, &shift, 16);
		unsigned int log2_max_mv_length_vertical = get_ue(CPB, &shift, 15);
		printf("<li>tiles_fixed_structure_flag: <code>%x</code></li>\n"
			"<li>motion_vectors_over_pic_boundaries_flag: <code>%x</code></li>\n"
			"<li>restricted_ref_pic_lists_flag: <code>%x</code></li>\n"
			"<li>min_spatial_segmentation_idc: <code>%u</code></li>\n"
			"<li>max_bytes_per_pic_denom: <code>%u</code></li>\n"
			"<li>max_bits_per_min_cu_denom: <code>%u</code></li>\n"
			"<li>log2_max_mv_length_horizontal: <code>%u</code></li>\n"
			"<li>log2_max_mv_length_vertical: <code>%u</code></li>\n",
			tiles_fixed_structure_flag,
			motion_vectors_over_pic_boundaries_flag,
			restricted_ref_pic_lists_flag,
			min_spatial_segmentation_idc,
			max_bytes_per_pic_denom,
			max_bits_per_min_cu_denom,
			log2_max_mv_length_horizontal,
			log2_max_mv_length_vertical);
	}
	return shift;
}



/**
 * Prints out the byte-aligned profile_tier_level() for the highest sub-layer,
 * and returns the pointer past it. Overflows for at most 86 set bytes.
 */
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
	c += 12;
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
	/* Sub-layer selection should occur at the demux level, hence any such info is ignored. */
	if (p->max_sub_layers > 1) {
		unsigned int sub_layer_flags = (c[0] << 8) | c[1];
		c += 2 + 11 * __builtin_popcount(sub_layer_flags & 0xaaa0) +
			__builtin_popcount(sub_layer_flags & 0x5550);
	}
	return c;
}



/**
 * Parses the SPS into a Rage265_parameter_set structure, and returns NULL.
 * Overflows for at most 2310 set bits.
 */
static const Rage265_picture *parse_sequence_parameter_set(Rage265_ctx *r, unsigned int lim) {
	static const char * const chroma_format_idc_names[4] = {"4:0:0", "4:2:0", "4:2:2", "4:4:4"};
	
	Rage265_parameter_set s = {0};
	s.num_tile_columns = 1;
	s.num_tile_rows = 1;
	s.loop_filter_across_tiles_enabled_flag = 1;
	unsigned int sps_video_parameter_set_id = r->CPB[0] >> 4;
	s.max_sub_layers = umin((r->CPB[0] >> 1) & 0x7, 6) + 1;
	s.temporal_id_nesting_flag = r->CPB[0] & 1;
	// shift reaches lim here in the worst case of overflow
	printf("<li>sps_video_parameter_set_id: <code>%u</code></li>\n"
		"<li>sps_max_sub_layers: <code>%u</code></li>\n"
		"<li>sps_temporal_id_nesting_flag: <code>%x</code></li>\n",
		sps_video_parameter_set_id,
		s.max_sub_layers,
		s.temporal_id_nesting_flag);
	unsigned int shift = 8 * (parse_profile_tier_level(&s, r->CPB + 1) - r->CPB);
	unsigned int sps_seq_parameter_set_id = get_ue(r->CPB, &shift, 15);
	s.ChromaArrayType = s.chroma_format_idc = get_ue(r->CPB, &shift, 3);
	printf("<li%s>sps_seq_parameter_set_id: <code>%u</code></li>\n"
		"<li>chroma_format_idc: <code>%u (%s)</code></li>\n",
		red_if(sps_seq_parameter_set_id > 0), sps_seq_parameter_set_id,
		s.chroma_format_idc, chroma_format_idc_names[s.chroma_format_idc]);
	if (s.chroma_format_idc == 3) {
		s.separate_colour_plane_flag = get_u1(r->CPB, &shift);
		s.ChromaArrayType &= s.separate_colour_plane_flag - 1;
		printf("<li>separate_colour_plane_flag: <code>%x</code></li>\n",
			s.separate_colour_plane_flag);
	}
	s.pic_width_in_luma_samples = umax(get_ue(r->CPB, &shift, 16888), 1);
	s.pic_height_in_luma_samples = umax(get_ue(r->CPB, &shift, 16888), 1);
	printf("<li>pic_width_in_luma_samples: <code>%u</code></li>\n"
		"<li>pic_height_in_luma_samples: <code>%u</code></li>\n",
		s.pic_width_in_luma_samples,
		s.pic_height_in_luma_samples);
	if (get_u1(r->CPB, &shift)) {
		unsigned int shiftX = (s.ChromaArrayType == 1 || s.ChromaArrayType == 2);
		unsigned int shiftY = (s.ChromaArrayType == 1);
		unsigned int limX = (s.pic_width_in_luma_samples - 1) >> shiftX;
		unsigned int limY = (s.pic_height_in_luma_samples - 1) >> shiftY;
		s.conf_win_left_offset = umin(get_raw_ue(r->CPB, &shift, 16887), limX) << shiftX;
		s.conf_win_right_offset = umin(get_raw_ue(r->CPB, &shift, 16887), limX - s.conf_win_left_offset) << shiftX;
		s.conf_win_top_offset = umin(get_raw_ue(r->CPB, &shift, 16887), limY) << shiftY;
		s.conf_win_bottom_offset = umin(get_raw_ue(r->CPB, &shift, 16887), limY - s.conf_win_top_offset) << shiftY;
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
	s.CtbLog2SizeY = umin(umax(s.MinCbLog2SizeY + get_raw_ue(r->CPB, &shift, 3), 4), 6);
	s.PicWidthInCtbsY = s.pic_width_in_luma_samples >> s.CtbLog2SizeY;
	s.PicHeightInCtbsY = s.pic_height_in_luma_samples >> s.CtbLog2SizeY;
	s.Log2MinTrafoSize = umin(get_raw_ue(r->CPB, &shift, 3) + 2,
		umin(s.MinCbLog2SizeY, 5));
	s.Log2MaxTrafoSize = umin(s.Log2MinTrafoSize + get_raw_ue(r->CPB, &shift, 3),
		umin(s.CtbLog2SizeY, 5));
	s.max_transform_hierarchy_depth_inter = umin(get_raw_ue(r->CPB, &shift, 4),
		s.CtbLog2SizeY - s.Log2MinTrafoSize);
	s.max_transform_hierarchy_depth_intra = umin(get_raw_ue(r->CPB, &shift, 4),
		s.CtbLog2SizeY - s.Log2MinTrafoSize);
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
		if (!get_u1(r->CPB, &shift)) {
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
		s.PcmBitDepth_Y = umin((u >> 4) + 1, s.BitDepth_Y);
		s.PcmBitDepth_C = umin((u & 0xf) + 1, s.BitDepth_C);
		s.Log2MinIpcmCbSizeY = umin(umax(get_raw_ue(r->CPB, &shift, 2) + 3,
			s.MinCbLog2SizeY), umin(s.CtbLog2SizeY, 5));
		s.Log2MaxIpcmCbSizeY = umin(s.Log2MinIpcmCbSizeY + get_raw_ue(r->CPB, &shift, 2),
			umin(s.CtbLog2SizeY, 5));
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
	printf("<li>sps_temporal_mvp_enabled_flag: <code>%x</code></li>\n"
		"<li>strong_intra_smoothing_enabled_flag: <code>%x</code></li>\n",
		s.temporal_mvp_enabled_flag,
		s.strong_intra_smoothing_enabled_flag);
	if (get_u1(r->CPB, &shift))
		shift = parse_vui_parameters(r->CPB, shift, s.max_sub_layers);
	unsigned int sps_extension_flag = get_u1(r->CPB, &shift);
	printf("<li>sps_extension_flag: <code>%x</code></li>\n", sps_extension_flag);
	if (sps_extension_flag && shift < lim)
		shift = lim;
	if (shift != lim)
		printf("<li style=\"color: red\">Bitstream overflow (%d bits)</li>\n", shift - lim);
	if (shift != lim || sps_seq_parameter_set_id > 0 || s.general_profile_space > 0)
		return NULL;
	
	/* Clear the r->CPB and reallocate the DPB when the image format changes. */
	s.pic_width_in_luma_samples &= -1 << s.MinCbLog2SizeY;
	s.pic_height_in_luma_samples &= -1 << s.MinCbLog2SizeY;
	s.image_offsets[1] = s.pic_width_in_luma_samples * s.pic_height_in_luma_samples;
	if (s.image_offsets[1] > 35651584) {
		s.pic_height_in_luma_samples = 35651584 / s.pic_width_in_luma_samples;
		s.image_offsets[1] = s.pic_width_in_luma_samples * s.pic_height_in_luma_samples;
	}
	unsigned int PicSizeInSamplesC = s.image_offsets[1] / 4 *
		(1 << s.chroma_format_idc >> 1);
	s.image_offsets[2] = s.image_offsets[1] + PicSizeInSamplesC;
	s.image_offsets[3] = s.image_offsets[2] + PicSizeInSamplesC;
	if ((s.chroma_format_idc ^ r->SPS.chroma_format_idc) |
		(s.pic_width_in_luma_samples ^ r->SPS.pic_width_in_luma_samples) |
		(s.pic_height_in_luma_samples ^ r->SPS.pic_height_in_luma_samples) |
		(s.BitDepth_Y ^ r->SPS.BitDepth_Y) | (s.BitDepth_C ^ r->SPS.BitDepth_C) |
		(s.max_dec_pic_buffering ^ r->SPS.max_dec_pic_buffering)) {
		free((uint8_t *)r->CPB);
		r->CPB = NULL;
		r->CPB_size = 0;
		size_t ctb_size = s.PicWidthInCtbsY * s.PicHeightInCtbsY * sizeof(Rage265_ctb);
		if (r->DPB[0].image != NULL)
			free(r->DPB[0].image);
		void *p = calloc(s.max_dec_pic_buffering, s.image_offsets[3] * 2 + ctb_size);
		for (unsigned int i = 0; i < s.max_dec_pic_buffering; i++) {
			r->DPB[i].image = p + i * (s.image_offsets[3] * 2 + ctb_size);
			r->DPB[i].CTBs = (Rage265_ctb *)(r->DPB[i].image + s.image_offsets[3]);
			r->DPB[i].needed_for_output = 0;
			r->DPB[i].grey_picture = 0;
			r->DPB[i].PicOrderCntVal = INT32_MIN;
		}
		memset(r->PPSs, 0, sizeof(r->PPSs));
	}
	r->SPS = s;
	memcpy(r->short_term_RPSs, short_term_RPSs, sizeof(short_term_RPSs));
	return NULL;
}



/** Prints out the VPS and returns NULL. Overflows for at most 801 set bits. */
static const Rage265_picture *parse_video_parameter_set(Rage265_ctx *r, unsigned int lim) {
	Rage265_parameter_set v;
	unsigned int u = htobe16(*(uint16_t *)r->CPB);
	unsigned int vps_video_parameter_set_id = u >> 12;
	unsigned int vps_max_layers = ((u >> 4) & 0x3f) + 1;
	v.max_sub_layers = umin((u >> 1) & 0x7, 6) + 1;
	unsigned int vps_temporal_id_nesting_flag = u & 1;
	// shift reaches lim here in the worst case of overflow
	printf("<li>vps_video_parameter_set_id: <code>%u</code></li>\n"
		"<li>vps_max_layers: <code>%u</code></li>\n"
		"<li>vps_max_sub_layers: <code>%u</code></li>\n"
		"<li>vps_temporal_id_nesting_flag: <code>%x</code></li>\n",
		vps_video_parameter_set_id,
		vps_max_layers,
		v.max_sub_layers,
		vps_temporal_id_nesting_flag);
	unsigned int shift = 8 * (parse_profile_tier_level(&v, r->CPB + 4) - r->CPB);
	unsigned int vps_sub_layer_ordering_info_present_flag = get_u1(r->CPB, &shift);
	for (unsigned int i = (v.max_sub_layers - 1) & -vps_sub_layer_ordering_info_present_flag; i < v.max_sub_layers; i++) {
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
	for (unsigned int i = 1; i < vps_num_layer_sets && shift < lim; i++) {
		printf("<li>layer_id_included_flags[%u]: <code>", i);
		for (unsigned int j = 0; j <= vps_max_layer_id; j++) {
			unsigned int layer_id_included_flag = get_u1(r->CPB, &shift);
			printf("%x", layer_id_included_flag);
		}
		printf("</code></li>\n");
	}
	if (get_u1(r->CPB, &shift)) {
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
		for (unsigned int i = 0; i < vps_num_hrd_parameters && shift < lim; i++) {
			unsigned int hrd_layer_set_idx = get_ue(r->CPB, &shift, 1023);
			printf("<ul>\n"
				"<li>hrd_layer_set_idx[%u]: <code>%u</code></li>\n",
				i, hrd_layer_set_idx);
			unsigned int cprms_present_flag = 0;
			if (i > 0)
				cprms_present_flag = get_u1(r->CPB, &shift);
			shift = parse_hrd_parameters(r->CPB, shift, cprms_present_flag, v.max_sub_layers);
			printf("</ul>\n");
		}
	}
	unsigned int vps_extension_flag = get_u1(r->CPB, &shift);
	printf("<li>vps_extension_flag: <code>%x</code></li>\n", vps_extension_flag);
	if (vps_extension_flag && shift < lim)
		shift = lim;
	if (shift != lim)
		printf("<li style=\"color: red\">Bitstream overflow (%d bits)</li>\n", shift - lim);
	return NULL;
}



/** Find the start of the next 00 00 0n pattern, returning len if none was found. */
#ifdef __SSSE3__
size_t Rage265_find_start_code(const uint8_t *buf, size_t len, unsigned int n) {
	ptrdiff_t chunk = (uint8_t *)((uintptr_t)buf & -sizeof(__m128i)) - buf;
	for (size_t u = 0; chunk < (ptrdiff_t)len; u = chunk += sizeof(__m128i)) {
		/* Skip chunks without a zero odd byte. */
		if ((_mm_movemask_epi8(_mm_cmpeq_epi8(*(__m128i *)(buf + chunk),
			_mm_setzero_si128())) & 0xaaaa) == 0)
			continue;
		size_t lim = umin(chunk + sizeof(__m128i) + 2, len);
		for (unsigned int start_code = -1; u < lim; u++) {
			start_code = ((start_code & 0xffff) << 8) | buf[u];
			if (start_code == n)
				return u - 2;
		}
	}
	return len;
}
#endif



/** Parses a NAL unit, returning the next picture to output or NULL. */
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
	typedef const Rage265_picture *(*Parser)(Rage265_ctx *, unsigned int);
	static const Parser parse_nal_unit[64] = {
		[0 ... 21] = parse_slice_segment_layer,
		[32] = parse_video_parameter_set,
		[33] = parse_sequence_parameter_set,
		[34] = parse_picture_parameter_set,
		[35] = parse_access_unit_delimiter,
		[36] = parse_end_of_seq,
		[37] = parse_end_of_bitstream,
	};
	
	/* Allocate the CPB. */
	if (len < 2)
		return NULL;
	const unsigned int suffix_size = 289; // largest overflow for a SPS
	size_t CPB_size = len - 2 + suffix_size;
	if (CPB_size > 800000000 / 8 + suffix_size) {
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
		dst--;
	unsigned int lim = 8 * (dst - r->CPB) + 7 - __builtin_ctz(*dst);
	memset(dst + 1, 0xff, suffix_size);
	
	/* Parse the nal_unit_header() and branch on nal_unit_type. */
	unsigned int nal_unit_header = (buf[0] << 8) | buf[1];
	r->nal_unit_type = nal_unit_header >> 9;
	unsigned int nuh_layer_id = (nal_unit_header >> 3) & 0x3f;
	r->TemporalId = (nal_unit_header - 1) & 0x7;
	printf("<ul class=\"frame\">\n"
		"<li%s>nal_unit_type: <code>%u (%s)</code></li>\n"
		"<li>nuh_layer_id: <code>%u</code></li>\n"
		"<li>TemporalId: <code>%u</code></li>\n",
		red_if(parse_nal_unit[r->nal_unit_type] == NULL), r->nal_unit_type, nal_unit_type_names[r->nal_unit_type],
		nuh_layer_id,
		r->TemporalId);
	const Rage265_picture *output = NULL;
	if (nuh_layer_id == 0 && parse_nal_unit[r->nal_unit_type] != NULL)
		output = parse_nal_unit[r->nal_unit_type](r, lim);
	printf("</ul>\n");
	return output;
}
