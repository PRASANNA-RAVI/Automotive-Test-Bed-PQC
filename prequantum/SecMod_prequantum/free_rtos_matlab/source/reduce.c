#include <stdint.h>
#include "params.h"
#include "reduce.h"

static const uint32_t qinv = 7679; // -inverse_mod(q,2^18)
static const uint32_t rlog = 18;

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a with 0 <= a <= Q*2^32,
*              compute r \equiv a*2^{-32} (mod Q) such that 0 < r < 2*Q.
*
* Arguments:   - uint64_t: element a
*
* Returns r.
**************************************************/
uint32_t montgomery_reduce(uint64_t a) {
  const uint64_t qinv = QINV;
  uint64_t t;

  t = a * qinv;
  t &= (1ULL << 32) - 1;
  t *= Q;
  t = a + t;
  return t >> 32;
}

/*************************************************
* Name:        reduce32
*
* Description: For finite field element a, compute r \equiv a (mod Q)
*              such that 0 <= r < 2*Q.
*
* Arguments:   - uint32_t: element a
*
* Returns r.
**************************************************/
uint32_t reduce32(uint32_t a) {
  uint32_t t;

  t = a & 0x7FFFFF;
  a >>= 23;
  t += ((a << 13) - a);
  return t;
}

/*************************************************
* Name:        freeze
*
* Description: For finite field element a, compute standard
*              representative r = a mod Q.
*
* Arguments:   - uint32_t: element a
*
* Returns r.
**************************************************/
uint32_t freeze(uint32_t a) {
  a = reduce32(a);
  a -= Q;
  a += ((int32_t)a >> 31) & Q;
  return a;
}

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q,
*              where R=2^18 (see value of rlog)
*
* Arguments:   - uint32_t a: input unsigned integer to be reduced; has to be in {0,...,2281446912}
*
* Returns:     unsigned integer in {0,...,2^13-1} congruent to a * R^-1 modulo q.
**************************************************/
uint16_t montgomery_reduce_kyber(uint32_t a)
{
  uint32_t u;

  u = (a * qinv);
  u &= ((1<<rlog)-1);
  u *= KYBER_Q;
  a = a + u;
  return a >> rlog;
}


/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q in {0,...,11768}
*
* Arguments:   - uint16_t a: input unsigned integer to be reduced
*
* Returns:     unsigned integer in {0,...,11768} congruent to a modulo q.
**************************************************/
uint16_t barrett_reduce(uint16_t a)
{
  uint32_t u;

  u = a >> 13;//((uint32_t) a * sinv) >> 16;
  u *= KYBER_Q;
  a -= u;
  return a;
}

/*************************************************
* Name:        freeze
*
* Description: Full reduction; given a 16-bit integer a, computes
*              unsigned integer a mod q.
*
* Arguments:   - uint16_t x: input unsigned integer to be reduced
*
* Returns:     unsigned integer in {0,...,q-1} congruent to a modulo q.
**************************************************/
uint16_t freeze_kyber(uint16_t x)
{
  uint16_t m,r;
  int16_t c;
  r = barrett_reduce(x);

  m = r - KYBER_Q;
  c = m;
  c >>= 15;
  r = m ^ ((r^m)&c);

  return r;
}
