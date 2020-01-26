#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"
#include "fips202.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  uint16_t coeffs[KYBER_N];
} poly_kyber;

void poly_compress(unsigned char *r, const poly_kyber *a);
void poly_decompress(poly_kyber *r, const unsigned char *a);

void poly_tobytes(unsigned char *r, const poly_kyber *a);
void poly_frombytes(poly_kyber *r, const unsigned char *a);

void poly_frommsg(poly_kyber *r, const unsigned char msg[KYBER_SYMBYTES]);
void poly_tomsg(unsigned char msg[KYBER_SYMBYTES], const poly_kyber *r);

void poly_getnoise(poly_kyber *r,const unsigned char *seed, unsigned char nonce);

void poly_ntt_kyber(poly_kyber *r);
void poly_invntt_kyber(poly_kyber *r);

void poly_add_kyber(poly_kyber *r, const poly_kyber *a, const poly_kyber *b);
void poly_sub_kyber(poly_kyber *r, const poly_kyber *a, const poly_kyber *b);

typedef struct {
  uint32_t coeffs[N];
} poly __attribute__((aligned(32)));

void poly_copy(poly *b, const poly *a);
void poly_freeze(poly *a);

void poly_add(poly *c, const poly *a, const poly *b);
void poly_sub(poly *c, const poly *a, const poly *b);
void poly_neg(poly *a);
void poly_shiftl(poly *a, unsigned int k);

void poly_ntt(poly *a);
void poly_invntt_montgomery(poly *a);
void poly_pointwise_invmontgomery(poly *c, const poly *a, const poly *b);

int  poly_chknorm(const poly *a, uint32_t B);
void poly_uniform(poly *a, unsigned char *buf);
void poly_uniform_eta(poly *a,
                      const unsigned char seed[SEEDBYTES],
                      unsigned char nonce);
void poly_uniform_gamma1m1(poly *a,
                           const unsigned char seed[SEEDBYTES + CRHBYTES],
                           uint16_t nonce);

void polyeta_pack(unsigned char *r, const poly *a);
void polyeta_unpack(poly *r, const unsigned char *a);

void polyt1_pack(unsigned char *r, const poly *a);
void polyt1_unpack(poly *r, const unsigned char *a);

void polyt0_pack(unsigned char *r, const poly *a);
void polyt0_unpack(poly *r, const unsigned char *a);

void polyz_pack(unsigned char *r, const poly *a);
void polyz_unpack(poly *r, const unsigned char *a);

void polyw1_pack(unsigned char *r, const poly *a);
#endif
