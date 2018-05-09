#ifndef _BAM_UTILS_H
#define _BAM_UTILS_H
#include "defs.h"
#include <cinttypes>

static inline void read_int32_t(int32_t& out, uint8_t* in, uint64_t& pos)
{
	out = 0;
	for (int j = 0; j < 4; ++j)
		out |= (uint32_t)in[pos++] << (j * 8);
}

static inline void read_uint32_t(uint32_t& out, uint8_t* in, uint64_t& pos)
{
	out = 0;
	for (int j = 0; j < 4; ++j)
		out |= (uint32_t)in[pos++] << (j * 8);
}

static inline void read_uint16_t(uint16_t& out, uint8_t* in, uint64_t& pos)
{
	out = 0;
	for (int j = 0; j < 2; ++j)
		out |= (uint16_t)in[pos++] << (j * 8);
}
#endif