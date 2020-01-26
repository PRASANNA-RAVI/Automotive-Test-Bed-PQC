#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

uint16_t freeze_kyber(uint16_t x);

uint16_t montgomery_reduce_kyber(uint32_t a);

uint16_t barrett_reduce(uint16_t a);

#define MONT 4193792U // 2^32 % Q
#define QINV 4236238847U // -q^(-1) mod 2^32

/* a <= Q*2^32 = > r < 2*Q */
uint32_t montgomery_reduce(uint64_t a);

/* r < 2*Q */
uint32_t reduce32(uint32_t a);

/* r < Q */
uint32_t freeze(uint32_t a);

#endif
