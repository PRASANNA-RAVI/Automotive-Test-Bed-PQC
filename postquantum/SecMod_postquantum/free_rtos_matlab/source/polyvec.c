#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "rounding.h"
#include "poly.h"
#include "polyvec.h"
#include "fips202.h"
#include "cbd.h"
#include "reduce.h"

#if (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 352))

/*************************************************
* Name:        polyvec_compress
*
* Description: Compress and serialize vector of polynomials
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_compress(unsigned char *r, const polyvec *a)
{
  int i,j,k;
  uint16_t t[8];
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N/8;j++)
    {
      for(k=0;k<8;k++)
        t[k] = ((((uint32_t)freeze_kyber(a->vec[i].coeffs[8*j+k]) << 11) + KYBER_Q/2)/ KYBER_Q) & 0x7ff;

      r[11*j+ 0] =  t[0] & 0xff;
      r[11*j+ 1] = (t[0] >>  8) | ((t[1] & 0x1f) << 3);
      r[11*j+ 2] = (t[1] >>  5) | ((t[2] & 0x03) << 6);
      r[11*j+ 3] = (t[2] >>  2) & 0xff;
      r[11*j+ 4] = (t[2] >> 10) | ((t[3] & 0x7f) << 1);
      r[11*j+ 5] = (t[3] >>  7) | ((t[4] & 0x0f) << 4);
      r[11*j+ 6] = (t[4] >>  4) | ((t[5] & 0x01) << 7);
      r[11*j+ 7] = (t[5] >>  1) & 0xff;
      r[11*j+ 8] = (t[5] >>  9) | ((t[6] & 0x3f) << 2);
      r[11*j+ 9] = (t[6] >>  6) | ((t[7] & 0x07) << 5);
      r[11*j+10] = (t[7] >>  3);
    }
    r += 352;
  }
}

/*************************************************
* Name:        polyvec_decompress
*
* Description: De-serialize and decompress vector of polynomials;
*              approximate inverse of polyvec_compress
*
* Arguments:   - polyvec *r:       pointer to output vector of polynomials
*              - unsigned char *a: pointer to input byte array
**************************************************/
void polyvec_decompress(polyvec *r, const unsigned char *a)
{
  int i,j;
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N/8;j++)
    {
      r->vec[i].coeffs[8*j+0] =  (((a[11*j+ 0]       | (((uint32_t)a[11*j+ 1] & 0x07) << 8)) * KYBER_Q) +1024) >> 11;
      r->vec[i].coeffs[8*j+1] = ((((a[11*j+ 1] >> 3) | (((uint32_t)a[11*j+ 2] & 0x3f) << 5)) * KYBER_Q) +1024) >> 11;
      r->vec[i].coeffs[8*j+2] = ((((a[11*j+ 2] >> 6) | (((uint32_t)a[11*j+ 3] & 0xff) << 2) |  (((uint32_t)a[11*j+ 4] & 0x01) << 10)) * KYBER_Q) + 1024) >> 11;
      r->vec[i].coeffs[8*j+3] = ((((a[11*j+ 4] >> 1) | (((uint32_t)a[11*j+ 5] & 0x0f) << 7)) * KYBER_Q) + 1024) >> 11;
      r->vec[i].coeffs[8*j+4] = ((((a[11*j+ 5] >> 4) | (((uint32_t)a[11*j+ 6] & 0x7f) << 4)) * KYBER_Q) + 1024) >> 11;
      r->vec[i].coeffs[8*j+5] = ((((a[11*j+ 6] >> 7) | (((uint32_t)a[11*j+ 7] & 0xff) << 1) |  (((uint32_t)a[11*j+ 8] & 0x03) <<  9)) * KYBER_Q) + 1024) >> 11;
      r->vec[i].coeffs[8*j+6] = ((((a[11*j+ 8] >> 2) | (((uint32_t)a[11*j+ 9] & 0x1f) << 6)) * KYBER_Q) + 1024) >> 11;
      r->vec[i].coeffs[8*j+7] = ((((a[11*j+ 9] >> 5) | (((uint32_t)a[11*j+10] & 0xff) << 3)) * KYBER_Q) + 1024) >> 11;
    }
    a += 352;
  }
}

#elif (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 320))

void polyvec_compress(unsigned char *r, const polyvec *a)
{
  int i,j,k;
  uint16_t t[4];
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N/4;j++)
    {
      for(k=0;k<4;k++)
        t[k] = ((((uint32_t)freeze_kyber(a->vec[i].coeffs[4*j+k]) << 10) + KYBER_Q/2)/ KYBER_Q) & 0x3ff;

      r[5*j+ 0] =  t[0] & 0xff;
      r[5*j+ 1] = (t[0] >>  8) | ((t[1] & 0x3f) << 2);
      r[5*j+ 2] = (t[1] >>  6) | ((t[2] & 0x0f) << 4);
      r[5*j+ 3] = (t[2] >>  4) | ((t[3] & 0x03) << 6);
      r[5*j+ 4] = (t[3] >>  2);
    }
    r += 320;
  }
}

void polyvec_decompress(polyvec *r, const unsigned char *a)
{
  int i,j;
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N/4;j++)
    {
      r->vec[i].coeffs[4*j+0] =  (((a[5*j+ 0]       | (((uint32_t)a[5*j+ 1] & 0x03) << 8)) * KYBER_Q) + 512) >> 10;
      r->vec[i].coeffs[4*j+1] = ((((a[5*j+ 1] >> 2) | (((uint32_t)a[5*j+ 2] & 0x0f) << 6)) * KYBER_Q) + 512) >> 10;
      r->vec[i].coeffs[4*j+2] = ((((a[5*j+ 2] >> 4) | (((uint32_t)a[5*j+ 3] & 0x3f) << 4)) * KYBER_Q) + 512) >> 10;
      r->vec[i].coeffs[4*j+3] = ((((a[5*j+ 3] >> 6) | (((uint32_t)a[5*j+ 4] & 0xff) << 2)) * KYBER_Q) + 512) >> 10;
    }
    a += 320;
  }
}

#elif (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 288))

void polyvec_compress(unsigned char *r, const polyvec *a)
{
  int i,j,k;
  uint16_t t[8];
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N/8;j++)
    {
      for(k=0;k<8;k++)
        t[k] = ((((uint32_t)freeze_kyber(a->vec[i].coeffs[8*j+k]) << 9) + KYBER_Q/2)/ KYBER_Q) & 0x1ff;

      r[9*j+ 0] =  t[0] & 0xff;
      r[9*j+ 1] = (t[0] >>  8) | ((t[1] & 0x7f) << 1);
      r[9*j+ 2] = (t[1] >>  7) | ((t[2] & 0x3f) << 2);
      r[9*j+ 3] = (t[2] >>  6) | ((t[3] & 0x1f) << 3);
      r[9*j+ 4] = (t[3] >>  5) | ((t[4] & 0x0f) << 4);
      r[9*j+ 5] = (t[4] >>  4) | ((t[5] & 0x07) << 5);
      r[9*j+ 6] = (t[5] >>  3) | ((t[6] & 0x03) << 6);
      r[9*j+ 7] = (t[6] >>  2) | ((t[7] & 0x01) << 7);
      r[9*j+ 8] = (t[7] >>  1);
    }
    r += 288;
  }
}

void polyvec_decompress(polyvec *r, const unsigned char *a)
{
  int i,j;
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N/8;j++)
    {
      r->vec[i].coeffs[8*j+0] =  (((a[9*j+ 0]       | (((uint32_t)a[9*j+ 1] & 0x01) << 8)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+1] = ((((a[9*j+ 1] >> 1) | (((uint32_t)a[9*j+ 2] & 0x03) << 7)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+2] = ((((a[9*j+ 2] >> 2) | (((uint32_t)a[9*j+ 3] & 0x07) << 6)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+3] = ((((a[9*j+ 3] >> 3) | (((uint32_t)a[9*j+ 4] & 0x0f) << 5)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+4] = ((((a[9*j+ 4] >> 4) | (((uint32_t)a[9*j+ 5] & 0x1f) << 4)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+5] = ((((a[9*j+ 5] >> 5) | (((uint32_t)a[9*j+ 6] & 0x3f) << 3)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+6] = ((((a[9*j+ 6] >> 6) | (((uint32_t)a[9*j+ 7] & 0x7f) << 2)) * KYBER_Q) + 256) >> 9;
      r->vec[i].coeffs[8*j+7] = ((((a[9*j+ 7] >> 7) | (((uint32_t)a[9*j+ 8] & 0xff) << 1)) * KYBER_Q) + 256) >> 9;
    }
    a += 288;
  }
}


#elif (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 256))

void polyvec_compress(unsigned char *r, const polyvec *a)
{
  int i,j,k;
  uint16_t t;
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N;j++)
    {
      r[j] = ((((uint32_t)freeze_kyber(a->vec[i].coeffs[j]) << 8) + KYBER_Q/2)/ KYBER_Q) & 0xff;
    }
    r += 256;
  }
}

void polyvec_decompress(polyvec *r, const unsigned char *a)
{
  int i,j;
  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_N;j++)
    {
      r->vec[i].coeffs[j] = ((a[j] * KYBER_Q) + 128) >> 8;
    }
    a += 256;
  }
}

#else
  #error "Unsupported compression of polyvec"
#endif

/*************************************************
* Name:        polyvec_tobytes
*
* Description: Serialize vector of polynomials
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_tobytes(unsigned char *r, const polyvec *a)
{
  int i;
  for(i=0;i<KYBER_K;i++)
    poly_tobytes(r+i*KYBER_POLYBYTES, &a->vec[i]);
}

/*************************************************
* Name:        polyvec_frombytes
*
* Description: De-serialize vector of polynomials;
*              inverse of polyvec_tobytes
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_frombytes(polyvec *r, const unsigned char *a)
{
  int i;
  for(i=0;i<KYBER_K;i++)
    poly_frombytes(&r->vec[i], a+i*KYBER_POLYBYTES);
}

/*************************************************
* Name:        polyvec_ntt
*
* Description: Apply forward NTT to all elements of a vector of polynomials
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_ntt(polyvec *r)
{
  int i;
  for(i=0;i<KYBER_K;i++)
    poly_ntt_kyber(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_invntt
*
* Description: Apply inverse NTT to all elements of a vector of polynomials
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_invntt(polyvec *r)
{
  int i;
  for(i=0;i<KYBER_K;i++)
    poly_invntt_kyber(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_pointwise_acc
*
* Description: Pointwise multiply elements of a and b and accumulate into r
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_pointwise_acc(poly_kyber *r, const polyvec *a, const polyvec *b)
{
  int i,j;
  uint16_t t;
  for(j=0;j<KYBER_N;j++)
  {
    t = montgomery_reduce_kyber(4613* (uint32_t)b->vec[0].coeffs[j]); // 4613 = 2^{2*18} % q
    r->coeffs[j] = montgomery_reduce_kyber(a->vec[0].coeffs[j] * t);
    for(i=1;i<KYBER_K;i++)
    {
      t = montgomery_reduce_kyber(4613* (uint32_t)b->vec[i].coeffs[j]);
      r->coeffs[j] += montgomery_reduce_kyber(a->vec[i].coeffs[j] * t);
    }
    r->coeffs[j] = barrett_reduce(r->coeffs[j]);
  }
}

/*************************************************
* Name:        polyvec_add
*
* Description: Add vectors of polynomials
*
* Arguments: - polyvec *r:       pointer to output vector of polynomials
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_add(polyvec *r, const polyvec *a, const polyvec *b)
{
  int i;
  for(i=0;i<KYBER_K;i++)
    poly_add_kyber(&r->vec[i], &a->vec[i], &b->vec[i]);

}

/**************************************************************/
/************ Vectors of polynomials of length L **************/
/**************************************************************/

/*************************************************
* Name:        polyvecl_freeze
*
* Description: Reduce coefficients of polynomials in vector of length L
*              to standard representatives.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
void polyvecl_freeze(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_freeze(v->vec+i);
}

/*************************************************
* Name:        polyvecl_add
*
* Description: Add vectors of polynomials of length L.
*              No modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - polyvecl *u: pointer to first summand
*              - polyvecl *v: pointer to second summand
**************************************************/
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_add(w->vec+i, u->vec+i, v->vec+i);
}

/*************************************************
* Name:        polyvecl_ntt
*
* Description: Forward NTT of all polynomials in vector of length L. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
void polyvecl_ntt(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_ntt(v->vec+i);
}

/*************************************************
* Name:        polyvecl_pointwise_acc_invmontgomery
*
* Description: Pointwise multiply vectors of polynomials of length L, multiply *              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*              Output coeffcient are less than 2*Q.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/
void polyvecl_pointwise_acc_invmontgomery(poly *w,
                                          const polyvecl *u,
                                          const polyvecl *v)
{
  unsigned int i;
  poly t;

  poly_pointwise_invmontgomery(w, u->vec+0, v->vec+0);

  for(i = 1; i < L; ++i) {
    poly_pointwise_invmontgomery(&t, u->vec+i, v->vec+i);
    poly_add(w, w, &t);
  }

  for(i = 0; i < N; ++i)
    w->coeffs[i] = reduce32(w->coeffs[i]);
}

/*************************************************
* Name:        polyvecl_chknorm
*
* Description: Check infinity norm of polynomials in vector of length L.
*              Assumes input coefficients to be standard representatives.
*
* Arguments:   - const polyvecl *v: pointer to vector
*              - uint32_t B: norm bound
*
* Returns 0 if norm of all polynomials is strictly smaller than B and 1
* otherwise.
**************************************************/
int polyvecl_chknorm(const polyvecl *v, uint32_t bound)  {
  unsigned int i;
  int ret = 0;

  for(i = 0; i < L; ++i)
    ret |= poly_chknorm(v->vec+i, bound);

  return ret;
}

/**************************************************************/
/************ Vectors of polynomials of length K **************/
/**************************************************************/

/*************************************************
* Name:        polyveck_freeze
*
* Description: Reduce coefficients of polynomials in vector of length K
*              to standard representatives.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_freeze(polyveck *v)  {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_freeze(v->vec+i);
}

/*************************************************
* Name:        polyveck_add
*
* Description: Add vectors of polynomials of length K.
*              No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - polyveck *u: pointer to first summand
*              - polyveck *v: pointer to second summand
**************************************************/
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_add(w->vec+i, u->vec+i, v->vec+i);
}

/*************************************************
* Name:        polyveck_sub
*
* Description: Subtract vectors of polynomials of length K.
*              Assumes coefficients of polynomials in input vectors to be less
*              than 2*Q. No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - polyveck *u: pointer to first input vector
*              - polyveck *v: pointer to second input vector to be subtracted
*                             from first input vector
**************************************************/
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_sub(w->vec+i, u->vec+i, v->vec+i);
}

/*************************************************
* Name:        polyveck_neg
*
* Description: Negate vector of polynomials of length K.
*              Assumes input coefficients to be less than 2*Q.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_neg(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_neg(v->vec+i);
}

/*************************************************
* Name:        polyveck_shiftl
*
* Description: Multiply vector of polynomials of Length K by 2^k.
*
* Arguments:   - polyveck *v: pointer to input/output vector
*              - unsigned int k: exponent
**************************************************/
void polyveck_shiftl(polyveck *v, unsigned int k) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_shiftl(v->vec+i, k);
}

/*************************************************
* Name:        polyveck_ntt
*
* Description: Forward NTT of all polynomials in vector of length K. Output *              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_ntt(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_ntt(v->vec+i);
}

/*************************************************
* Name:        polyveck_invntt_montgomery
*
* Description: Inverse NTT and multiplication by 2^{32} of polynomials
*              in vector of length K. Input coefficients need to be less
*              than 2*Q.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_invntt_montgomery(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_invntt_montgomery(v->vec+i);
}

/*************************************************
* Name:        polyveck_chknorm
*
* Description: Check infinity norm of polynomials in vector of length K.
*              Assumes input coefficients to be standard representatives.
*
* Arguments:   - const polyveck *v: pointer to vector
*              - uint32_t B: norm bound
*
* Returns 0 if norm of all polynomials is strictly smaller than B and 1
* otherwise.
**************************************************/
int polyveck_chknorm(const polyveck *v, uint32_t bound) {
  unsigned int i;
  int ret = 0;

  for(i = 0; i < K; ++i)
    ret |= poly_chknorm(v->vec+i, bound);

  return ret;
}

/*************************************************
* Name:        polyveck_power2round
*
* Description: For all coefficients a of polynomials in vector of length K,
*              compute a0, a1 such that a = a1*2^D + a0
*              with -2^{D/2} < a0 <= 2^{D/2}. Assumes a to be standard
*              representative.
*
* Arguments:   - polyveck *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyveck *v0: pointer to output vector of polynomials with
*                              coefficients a0
*              - polyveck *v: pointer to input vector
**************************************************/
void polyveck_power2round(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i, j;

  for(i = 0; i < K; ++i)
    for(j = 0; j < N; ++j)
      v1->vec[i].coeffs[j] = power2round(v->vec[i].coeffs[j],
                                         &v0->vec[i].coeffs[j]);
}

/*************************************************
* Name:        polyveck_decompose
*
* Description: For all coefficients a of polynomials in vector of length K,
*              compute high and low bits a0, a1 such a = a1*ALPHA + a0
*              with -ALPHA/2 < a0 <= ALPHA/2 except if a = Q-1 where
*              a1 = 0 and a0 = -1. Assumes a to be standard representative.
*
* Arguments:   - polyveck *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyveck *v0: pointer to output vector of polynomials with
*                              coefficients a0
*              - polyveck *v: pointer to input vector
**************************************************/
void polyveck_decompose(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i, j;

  for(i = 0; i < K; ++i)
    for(j = 0; j < N; ++j)
      v1->vec[i].coeffs[j] = decompose(v->vec[i].coeffs[j],
                                       &v0->vec[i].coeffs[j]);
}

/*************************************************
* Name:        polyveck_make_hint
*
* Description: Compute hint vector. The coefficients of the polynomials indicate
*              whether or not the high bits of the corresponding coefficients
*              of the input polynomials differ.
*
* Arguments:   - polyveck *h: pointer to output vector
*              - const polyveck *u: pointer to first input vector
*              - const polyveck *u: pointer to second input vector
*
* Returns number of 1 bits.
**************************************************/
unsigned int polyveck_make_hint(polyveck *h,
                                const polyveck *u,
                                const polyveck *v)
{
  unsigned int i, j, s = 0;

  for(i = 0; i < K; ++i)
    for(j = 0; j < N; ++j) {
      h->vec[i].coeffs[j] = make_hint(u->vec[i].coeffs[j], v->vec[i].coeffs[j]);
      s += h->vec[i].coeffs[j];
    }

  return s;
}

/*************************************************
* Name:        polyveck_use_hint
*
* Description: Use hint vector to correct the high bits of input vector.
*
* Arguments:   - polyveck *w: pointer to output vector of polynomials with
*                             corrected high bits
*              - polyveck *u: pointer to input vector
*              - polyveck *h: pointer to input hint vector
**************************************************/
void polyveck_use_hint(polyveck *w, const polyveck *u, const polyveck *h) {
  unsigned int i, j;

  for(i = 0; i < K; ++i)
    for(j = 0; j < N; ++j)
      w->vec[i].coeffs[j] = use_hint(u->vec[i].coeffs[j], h->vec[i].coeffs[j]);
}
