#include <stdint.h>
#include <stdio.h>
#include "api.h"
#include "params.h"
#include "sign.h"
// #include "randombytes.h"
#include "fips202.h"
#include "poly.h"
#include "polyvec.h"
#include "packing.h"

static void pack_pk(unsigned char pk[PK_SIZE_PACKED],
             const unsigned char rho[SEEDBYTES],
             const polyveck *t1)
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    pk[i] = rho[i];
  pk += SEEDBYTES;

  for(i = 0; i < K; ++i)
    polyt1_pack(pk + i*POLT1_SIZE_PACKED, t1->vec+i);
}

/*************************************************
* Name:        unpack_pk
*
* Description: Unpack public key pk = (rho, t1).
*
* Arguments:   - const unsigned char rho[]: output byte array for rho
*              - const polyveck *t1: pointer to output vector t1
*              - unsigned char pk[]: byte array containing bit-packed pk
**************************************************/
static void unpack_pk(unsigned char rho[SEEDBYTES],
               polyveck *t1,
               const unsigned char pk[PK_SIZE_PACKED])
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    rho[i] = pk[i];
  pk += SEEDBYTES;

  for(i = 0; i < K; ++i)
    polyt1_unpack(t1->vec+i, pk + i*POLT1_SIZE_PACKED);
}

/*************************************************
* Name:        pack_sk
*
* Description: Bit-pack secret key sk = (rho, key, tr, s1, s2, t0).
*
* Arguments:   - unsigned char sk[]: output byte array
*              - const unsigned char rho[]: byte array containing rho
*              - const unsigned char key[]: byte array containing key
*              - const unsigned char tr[]: byte array containing tr
*              - const polyvecl *s1: pointer to vector s1
*              - const polyveck *s2: pointer to vector s2
*              - const polyveck *t0: pointer to vector t0
**************************************************/
static void pack_sk(unsigned char sk[SK_SIZE_PACKED],
             const unsigned char rho[SEEDBYTES],
             const unsigned char key[SEEDBYTES],
             const unsigned char tr[CRHBYTES],
             const polyvecl *s1,
             const polyveck *s2,
             const polyveck *t0)
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    sk[i] = rho[i];
  sk += SEEDBYTES;

  for(i = 0; i < SEEDBYTES; ++i)
    sk[i] = key[i];
  sk += SEEDBYTES;

  for(i = 0; i < CRHBYTES; ++i)
    sk[i] = tr[i];
  sk += CRHBYTES;

  for(i = 0; i < L; ++i)
    polyeta_pack(sk + i*POLETA_SIZE_PACKED, s1->vec+i);
  sk += L*POLETA_SIZE_PACKED;

  for(i = 0; i < K; ++i)
    polyeta_pack(sk + i*POLETA_SIZE_PACKED, s2->vec+i);
  sk += K*POLETA_SIZE_PACKED;

  for(i = 0; i < K; ++i)
    polyt0_pack(sk + i*POLT0_SIZE_PACKED, t0->vec+i);
}

/*************************************************
* Name:        unpack_sk
*
* Description: Unpack secret key sk = (rho, key, tr, s1, s2, t0).
*
* Arguments:   - const unsigned char rho[]: output byte array for rho
*              - const unsigned char key[]: output byte array for key
*              - const unsigned char tr[]: output byte array for tr
*              - const polyvecl *s1: pointer to output vector s1
*              - const polyveck *s2: pointer to output vector s2
*              - const polyveck *r0: pointer to output vector t0
*              - unsigned char sk[]: byte array containing bit-packed sk
**************************************************/
static void unpack_sk(unsigned char rho[SEEDBYTES],
               unsigned char key[SEEDBYTES],
               unsigned char tr[CRHBYTES],
               polyvecl *s1,
               polyveck *s2,
               polyveck *t0,
               const unsigned char sk[SK_SIZE_PACKED])
{
  unsigned int i;

  for(i = 0; i < SEEDBYTES; ++i)
    rho[i] = sk[i];
  sk += SEEDBYTES;

  for(i = 0; i < SEEDBYTES; ++i)
    key[i] = sk[i];
  sk += SEEDBYTES;

  for(i = 0; i < CRHBYTES; ++i)
    tr[i] = sk[i];
  sk += CRHBYTES;

  for(i=0; i < L; ++i)
    polyeta_unpack(s1->vec+i, sk + i*POLETA_SIZE_PACKED);
  sk += L*POLETA_SIZE_PACKED;

  for(i=0; i < K; ++i)
    polyeta_unpack(s2->vec+i, sk + i*POLETA_SIZE_PACKED);
  sk += K*POLETA_SIZE_PACKED;

  for(i=0; i < K; ++i)
    polyt0_unpack(t0->vec+i, sk + i*POLT0_SIZE_PACKED);
}

/*************************************************
* Name:        pack_sig
*
* Description: Bit-pack signature sig = (z, h, c).
*
* Arguments:   - unsigned char sig[]: output byte array
*              - const polyvecl *z: pointer to vector z
*              - const polyveck *h: pointer to hint vector h
*              - const poly *c: pointer to challenge polynomial
**************************************************/
static void pack_sig(unsigned char sig[SIG_SIZE_PACKED],
              const polyvecl *z,
              const polyveck *h,
              const poly *c)
{
  unsigned int i, j, k;
  uint64_t signs, mask;

  for(i = 0; i < L; ++i)
    polyz_pack(sig + i*POLZ_SIZE_PACKED, z->vec+i);
  sig += L*POLZ_SIZE_PACKED;

  /* Encode h */
  k = 0;
  for(i = 0; i < K; ++i) {
    for(j = 0; j < N; ++j)
      if(h->vec[i].coeffs[j] == 1)
        sig[k++] = j;

    sig[OMEGA + i] = k;
  }
  while(k < OMEGA) sig[k++] = 0;
  sig += OMEGA + K;

  /* Encode c */
  signs = 0;
  mask = 1;
  for(i = 0; i < N/8; ++i) {
    sig[i] = 0;
    for(j = 0; j < 8; ++j) {
      if(c->coeffs[8*i+j] != 0) {
        sig[i] |= (1 << j);
        if(c->coeffs[8*i+j] == (Q - 1)) signs |= mask;
        mask <<= 1;
      }
    }
  }
  sig += N/8;
  for(i = 0; i < 8; ++i)
    sig[i] = signs >> 8*i;
}

/*************************************************
* Name:        unpack_sig
*
* Description: Unpack signature sig = (z, h, c).
*
* Arguments:   - const polyvecl *z: pointer to output vector z
*              - const polyveck *h: pointer to output hint vector h
*              - const poly *c: pointer to output challenge polynomial
*              - unsigned char sig[]: byte array containing bit-packed signature
**************************************************/
static void unpack_sig(polyvecl *z,
                polyveck *h,
                poly *c,
                const unsigned char sig[SIG_SIZE_PACKED])
{
  unsigned int i, j, k;
  uint64_t signs, mask;

  for(i = 0; i < L; ++i)
    polyz_unpack(z->vec+i, sig + i*POLZ_SIZE_PACKED);
  sig += L*POLZ_SIZE_PACKED;

  /* Decode h */
  k = 0;
  for(i = 0; i < K; ++i) {
    for(j = 0; j < N; ++j)
      h->vec[i].coeffs[j] = 0;

    for(j = k; j < sig[OMEGA + i]; ++j)
      h->vec[i].coeffs[sig[j]] = 1;

    k = sig[OMEGA + i];
  }
  sig += OMEGA + K;

  /* Decode c */
  for(i = 0; i < N; ++i)
    c->coeffs[i] = 0;

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)sig[N/8+i] << 8*i;

  mask = 1;
  for(i = 0; i < N/8; ++i) {
    for(j = 0; j < 8; ++j) {
      if((sig[i] >> j) & 0x01) {
        c->coeffs[8*i+j] = (signs & mask) ? Q - 1 : 1;
        mask <<= 1;
      }
    }
  }
}

static int shift_lfsr(unsigned int *lfsr, unsigned int polynomial_mask)
{
    int feedback;

    feedback = *lfsr & 1;
    *lfsr >>= 1;
    if(feedback == 1)
        *lfsr ^= polynomial_mask;
    return *lfsr;
}

static int get_random(void)
{
    int temp;
    unsigned int POLY_MASK_HERE_1 = 0xB765879A;
    unsigned int POLY_MASK_HERE_2 = 0x55BBEEFF;
    static unsigned int lfsr_1 = 0x55AAEEFF;
    static unsigned int lfsr_2 = 0xFFAA8844;
    shift_lfsr(&lfsr_1, POLY_MASK_HERE_1);
    shift_lfsr(&lfsr_2, POLY_MASK_HERE_2);
    temp = (shift_lfsr(&lfsr_1, POLY_MASK_HERE_1) ^ shift_lfsr(&lfsr_2, POLY_MASK_HERE_2)) & 0XFF;
    // printf("%02x\n",temp);
    return (temp);
}

static void get_random_bits(unsigned char *x, int length)
{
    unsigned int i;
    // for(i=0;i<length;i++)
    // {
    //   if(i == 0)
    //       rand_value = get_random();
    //   else
    //       rand_value = (rand_value<<8) | get_random();
    // }
    for(i=0;i<length;i++)
    {
        *(x+i) = get_random();
    }
}

/*************************************************
* Name:        expand_mat
*
* Description: Implementation of ExpandA. Generates matrix A with uniformly
*              random coefficients a_{i,j} by performing rejection
*              sampling on the output stream of SHAKE128(rho|i|j).
*
* Arguments:   - polyvecl mat[K]: output matrix
*              - const unsigned char rho[]: byte array containing seed rho
**************************************************/
void expand_mat(polyvecl mat[K], const unsigned char rho[SEEDBYTES]) {
  unsigned int i, j, pos, ctr;
  unsigned char inbuf[SEEDBYTES + 1];
  /* Don't change this to smaller values,
   * sampling later assumes sufficient SHAKE output!
   * Probability that we need more than 5 blocks: < 2^{-132}.
   * Probability that we need more than 6 blocks: < 2^{-546}. */
  unsigned char outbuf[5*SHAKE128_RATE];
  uint32_t val;

  for(i = 0; i < SEEDBYTES; ++i)
    inbuf[i] = rho[i];

  for(i = 0; i < K; ++i) {
    for(j = 0; j < L; ++j) {
      ctr = pos = 0;
      inbuf[SEEDBYTES] = i + (j << 4);

      shake128(outbuf, sizeof(outbuf), inbuf, SEEDBYTES + 1);

      while(ctr < N) {
        val  = outbuf[pos++];
        val |= (uint32_t)outbuf[pos++] << 8;
        val |= (uint32_t)outbuf[pos++] << 16;
        val &= 0x7FFFFF;

        /* Rejection sampling */
        if(val < Q)
          mat[i].vec[j].coeffs[ctr++] = val;
      }
    }
  }
}

/*************************************************
* Name:        challenge
*
* Description: Implementation of H. Samples polynomial with 60 nonzero
*              coefficients in {-1,1} using the output stream of
*              SHAKE256(mu|w1).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const unsigned char mu[]: byte array containing mu
*              - const polyveck *w1: pointer to vector w1
**************************************************/
void challenge(poly *c,
               const unsigned char mu[CRHBYTES],
               const polyveck *w1)
{
  unsigned int i, b, pos;
  unsigned char inbuf[CRHBYTES + K*POLW1_SIZE_PACKED];
  unsigned char outbuf[SHAKE256_RATE];
  uint64_t state[25], signs, mask;

  for(i = 0; i < CRHBYTES; ++i)
    inbuf[i] = mu[i];
  for(i = 0; i < K; ++i)
    polyw1_pack(inbuf + CRHBYTES + i*POLW1_SIZE_PACKED, w1->vec+i);

  shake256_absorb(state, inbuf, sizeof(inbuf));
  shake256_squeezeblocks(outbuf, 1, state);

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)outbuf[i] << 8*i;

  pos = 8;
  mask = 1;

  for(i = 0; i < N; ++i)
    c->coeffs[i] = 0;

  for(i = 196; i < 256; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(outbuf, 1, state);
        pos = 0;
      }

      b = outbuf[pos++];
    } while(b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = (signs & mask) ? Q - 1 : 1;
    mask <<= 1;
  }
}

/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - unsigned char *pk: pointer to output public key (allocated
*                                   array of CRYPTO_PUBLICKEYBYTES bytes)
*              - unsigned char *sk: pointer to output private key (allocated
*                                   array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
  unsigned int i;
  unsigned char seedbuf[3*SEEDBYTES];
  unsigned char tr[CRHBYTES];
  unsigned char *rho, *rhoprime, *key;
  uint16_t nonce = 0;
  polyvecl mat[K];
  polyvecl s1, s1hat;
  polyveck s2, t, t1, t0;

  /* Expand 32 bytes of randomness into rho, rhoprime and key */
  // randombytes(seedbuf, SEEDBYTES);
  get_random_bits(seedbuf, SEEDBYTES);
  // for(i=0;i<SEEDBYTES;i++)
  //   seedbuf[i] = 0xAA;
  shake256(seedbuf, 3*SEEDBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime = rho + SEEDBYTES;
  key = rho + 2*SEEDBYTES;

  /* Expand matrix */
  expand_mat(mat, rho);

  /* Sample short vectors s1 and s2 */
  for(i = 0; i < L; ++i)
    poly_uniform_eta(&s1.vec[i], rhoprime, nonce++);
  for(i = 0; i < K; ++i)
    poly_uniform_eta(&s2.vec[i], rhoprime, nonce++);

  /* Matrix-vector multiplication */
  s1hat = s1;
  polyvecl_ntt(&s1hat);
  for(i = 0; i < K; ++i) {
    polyvecl_pointwise_acc_invmontgomery(&t.vec[i], mat+i, &s1hat);
    poly_invntt_montgomery(t.vec+i);
  }

  /* Add noise vector s2 */
  polyveck_add(&t, &t, &s2);

  /* Extract t1 and write public key */
  polyveck_freeze(&t);
  polyveck_power2round(&t1, &t0, &t);
  pack_pk(pk, rho, &t1);

  /* Compute CRH(rho, t1) and write secret key */
  shake256(tr, CRHBYTES, pk, CRYPTO_PUBLICKEYBYTES);
  pack_sk(sk, rho, key, tr, &s1, &s2, &t0);

  return 0;
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - unsigned char *sm: pointer to output signed message (allocated
*                                   array with CRYPTO_BYTES + mlen bytes)
*              - unsigned long long *smlen: pointer to output length of signed
*                                           message
*              - const unsigned char *m: pointer to message to be signed
*              - unsigned long long mlen: length of message
*              - const unsigned char *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign(unsigned char *sm,
                unsigned long long *smlen,
                const unsigned char *m,
                unsigned long long mlen,
                const unsigned char *sk)
{
  unsigned long long i, j;
  unsigned int n;
  unsigned char seedbuf[2*SEEDBYTES + CRHBYTES];
  unsigned char *rho, *key, *mu, *tr;
  uint16_t nonce = 0;
  poly     c, chat;
  polyvecl mat[K], s1, y, yhat, z;
  polyveck s2, t0, w, w1;
  polyveck h, wcs2, wcs20, ct0, tmp;

  rho = seedbuf;
  key = seedbuf + SEEDBYTES;
  mu = seedbuf + 2*SEEDBYTES;
  tr = sm + CRYPTO_BYTES - CRHBYTES;
  unpack_sk(rho, key, tr, &s1, &s2, &t0, sk);

  /* Copy message at the end of the sm buffer */
  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + i] = m[i];

  /* Compute CRH(tr, msg) */
  shake256(mu, CRHBYTES, sm + CRYPTO_BYTES - CRHBYTES, CRHBYTES + mlen);

  /* Expand matrix and transform vectors */
  expand_mat(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);

  rej:
  /* Sample intermediate vector y */
  for(i = 0; i < L; ++i)
    poly_uniform_gamma1m1(y.vec+i, key, nonce++);

  /* Matrix-vector multiplication */
  yhat = y;
  polyvecl_ntt(&yhat);
  for(i = 0; i < K; ++i) {
    polyvecl_pointwise_acc_invmontgomery(w.vec+i, mat+i, &yhat);
    poly_invntt_montgomery(w.vec+i);
  }

  /* Decompose w and call the random oracle */
  polyveck_freeze(&w);
  polyveck_decompose(&w1, &tmp, &w);
  challenge(&c, mu, &w1);

  /* Compute z, reject if it reveals secret */
  chat = c;
  poly_ntt(&chat);
  for(i = 0; i < L; ++i) {
    poly_pointwise_invmontgomery(z.vec+i, &chat, s1.vec+i);
    poly_invntt_montgomery(z.vec+i);
  }
  polyvecl_add(&z, &z, &y);
  polyvecl_freeze(&z);
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    goto rej;

  /* Compute w - cs2, reject if w1 can not be computed from it */
  for(i = 0; i < K; ++i) {
    poly_pointwise_invmontgomery(wcs2.vec+i, &chat, s2.vec+i);
    poly_invntt_montgomery(wcs2.vec+i);
  }
  polyveck_sub(&wcs2, &w, &wcs2);
  polyveck_freeze(&wcs2);
  polyveck_decompose(&tmp, &wcs20, &wcs2);
  polyveck_freeze(&wcs20);
  if(polyveck_chknorm(&wcs20, GAMMA2 - BETA))
    goto rej;

  for(i = 0; i < K; ++i)
    for(j = 0; j < N; ++j)
      if(tmp.vec[i].coeffs[j] != w1.vec[i].coeffs[j])
        goto rej;

  /* Compute hints for w1 */
  for(i = 0; i < K; ++i) {
    poly_pointwise_invmontgomery(ct0.vec+i, &chat, t0.vec+i);
    poly_invntt_montgomery(ct0.vec+i);
  }

  polyveck_freeze(&ct0);
  if(polyveck_chknorm(&ct0, GAMMA2))
    goto rej;

  polyveck_add(&tmp, &wcs2, &ct0);
  polyveck_neg(&ct0);
  polyveck_freeze(&tmp);
  n = polyveck_make_hint(&h, &tmp, &ct0);
  if(n > OMEGA)
    goto rej;

  /* Write signature */
  pack_sig(sm, &z, &h, &c);

  *smlen = mlen + CRYPTO_BYTES;
  return 0;
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - unsigned char *m: pointer to output message (allocated
*                                  array with smlen bytes), can be equal to sm
*              - unsigned long long *mlen: pointer to output length of message
*              - const unsigned char *sm: pointer to signed message
*              - unsigned long long smlen: length of signed message
*              - const unsigned char *sk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open(unsigned char *m,
                     unsigned long long *mlen,
                     const unsigned char *sm,
                     unsigned long long smlen,
                     const unsigned char *pk)
{
  unsigned long long i;
  unsigned char rho[SEEDBYTES];
  unsigned char mu[CRHBYTES];
  poly     c, chat, cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h, tmp1, tmp2;

  if(smlen < CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES;

  unpack_pk(rho, &t1, pk);
  unpack_sig(&z, &h, &c, sm);
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    goto badsig;

  /* Compute CRH(CRH(rho, t1), msg) using m as "playground" buffer */
  for(i = 0; i < CRYPTO_PUBLICKEYBYTES; ++i)
    m[CRYPTO_BYTES - CRYPTO_PUBLICKEYBYTES + i] = pk[i];

  if(sm != m)
  {
    for(i = 0; i < *mlen; ++i)
    {
      m[CRYPTO_BYTES + i] = sm[CRYPTO_BYTES + i];
    }
  }
    // printf("\n");

  shake256(m + CRYPTO_BYTES - CRHBYTES, CRHBYTES,
           m + CRYPTO_BYTES - CRYPTO_PUBLICKEYBYTES, CRYPTO_PUBLICKEYBYTES);
  shake256(mu, CRHBYTES, m + CRYPTO_BYTES - CRHBYTES, CRHBYTES + *mlen);

  expand_mat(mat, rho);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  polyvecl_ntt(&z);
  for(i = 0; i < K ; ++i)
    polyvecl_pointwise_acc_invmontgomery(tmp1.vec+i, mat+i, &z);

  chat = c;
  poly_ntt(&chat);
  polyveck_shiftl(&t1, D);
  polyveck_ntt(&t1);
  for(i = 0; i < K; ++i)
    poly_pointwise_invmontgomery(tmp2.vec+i, &chat, t1.vec+i);

  polyveck_sub(&tmp1, &tmp1, &tmp2);
  polyveck_freeze(&tmp1);  // reduce32 would be sufficient
  polyveck_invntt_montgomery(&tmp1);

  /* Reconstruct w1 */
  polyveck_freeze(&tmp1);
  polyveck_use_hint(&w1, &tmp1, &h);

  /* Call random oracle and verify challenge */
  challenge(&cp, mu, &w1);
  for(i = 0; i < N; ++i)
    if(c.coeffs[i] != cp.coeffs[i])
      goto badsig;

  /* All good, copy msg, return 0 */
  for(i = 0; i < *mlen; ++i)
    m[i] = sm[CRYPTO_BYTES + i];

  return 0;

  /* Signature verification failed */
  badsig:
  *mlen = (unsigned long long) -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}
