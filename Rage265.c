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



static void parse_profile_tier_level(Rage265_worker *w) {
	
}



static void parse_hrd_parameters(Rage265_worker *w,
	unsigned int cprms_present_flag, unsigned int max_sub_layers) {
	
}



/**
 * This function parses the SPS into a Rage265_SPS structure, and if no error
 * was detected stores it into the Rage265_ctx object.
 */



/**
 * This function currently only prints the VPS to stdout.
 */
static unsigned int parse_VPS(Rage265_ctx *r, Rage265_worker *w) {
	unsigned int u = htobe16(*(uint16_t *)w->c.CPB);
	unsigned int vps_video_parameter_set_id = u >> 12;
	unsigned int vps_max_layers = ((u >> 4) & 0x3f) + 1;
	unsigned int vps_max_sub_layers = ((u >> 1) + 1) & 0x7;
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
			"<li>vps_max_latency_increase[%u]: <code>%d</code></li>\n"
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
	};
	
	/* On first call, initialise the main structure. */
	if (r->max_workers == 0)
		r->max_workers = get_nprocs();
	if (r->workers == NULL) {
		r->workers = calloc(r->max_workers, sizeof(Rage265_worker));
		if (r->workers == NULL)
			return 1 << RAGE265_ERROR_NO_MEMORY;
		r->lock = (pthread_mutex_t)PTHREAD_MUTEX_INITIALIZER;
		r->worker_available = (pthread_cond_t)PTHREAD_COND_INITIALIZER;
	}
	
	/* Assign an available worker to this NAL payload. */
	pthread_mutex_lock(&r->lock);
	Rage265_worker *w = r->workers;
	for (Rage265_worker *end = w += r->max_workers; w == end; ) {
		for (w = r->workers; w < end && w->target != NULL; w++);
		if (w == end)
			pthread_cond_wait(&r->worker_available, &r->lock);
	}
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
