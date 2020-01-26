#include <string.h>
#include "indcpa.h"
#include "poly.h"
#include "polyvec.h"
// #include "randombytes.h"
#include "fips202.h"
#include "ntt.h"
#include "packing.h"

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
    static unsigned int lfsr_1 = 0xAABBCCDD;
    static unsigned int lfsr_2 = 0x00AA5566;
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

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              SHAKE-128
*
* Arguments:   - polyvec *a:                pointer to ouptput matrix A
*              - const unsigned char *seed: pointer to input seed
*              - int transposed:            boolean deciding whether A or A^T is generated
**************************************************/
void gen_matrix(polyvec *a, const unsigned char *seed, int transposed) // Not static for benchmarking
{
  unsigned int pos=0, ctr;
  uint16_t val;
  unsigned int nblocks=4;
  int i,j;
  uint64_t state[25]; // SHAKE state
  unsigned char extseed[KYBER_SYMBYTES+2];
  unsigned char buff[672];

  for(i=0;i<KYBER_SYMBYTES;i++)
    extseed[i] = seed[i];


  for(i=0;i<KYBER_K;i++)
  {
    for(j=0;j<KYBER_K;j++)
    {
      ctr = pos = 0;
      if(transposed)
      {
        extseed[KYBER_SYMBYTES]   = i;
        extseed[KYBER_SYMBYTES+1] = j;
      }
      else
      {
        extseed[KYBER_SYMBYTES]   = j;
        extseed[KYBER_SYMBYTES+1] = i;
      }

      shake128_absorb(state,extseed,KYBER_SYMBYTES+2);
      shake128_squeezeblocks(buff,nblocks,state);

      while(ctr < KYBER_N)
      {
        val = (buff[pos] | ((uint16_t) buff[pos+1] << 8)) & 0x1fff;
        if(val < KYBER_Q)
        {
            a[i].vec[j].coeffs[ctr++] = val;
        }
        pos += 2;

        if(pos > SHAKE128_RATE*nblocks-2)
        {
          nblocks = 1;
          shake128_squeezeblocks(buff,nblocks,state);
          pos = 0;
        }
      }
    }
  }
}


/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - unsigned char *pk: pointer to output public key
*              - unsigned char *sk: pointer to output private key
**************************************************/
void indcpa_keypair(unsigned char *pk,
                   unsigned char *sk)
{
  polyvec a[KYBER_K], e, pkpv, skpv;
  unsigned char buf[KYBER_SYMBYTES+KYBER_SYMBYTES];
  unsigned char *publicseed = buf;
  unsigned char *noiseseed = buf+KYBER_SYMBYTES;
  int i;
  unsigned char nonce=0;

  // randombytes(buf, KYBER_SYMBYTES);
  get_random_bits(buf, KYBER_SYMBYTES);
  // for(i=0;i<KYBER_SYMBYTES;i++)
  //   buf[i] = 0xAA;
  sha3_512(buf, buf, KYBER_SYMBYTES);

  gen_a(a, publicseed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise(skpv.vec+i,noiseseed,nonce++);

    // printf("Public key\n");
    // for(i=0;i<KYBER_N;i++)
    // {
    //     printf("%d, ",skpv.vec[1].coeffs[i]);
    // }
    // printf("\n");

  polyvec_ntt(&skpv);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise(e.vec+i,noiseseed,nonce++);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_pointwise_acc(&pkpv.vec[i],&skpv,a+i);

  polyvec_invntt(&pkpv);
  polyvec_add(&pkpv,&pkpv,&e);

  pack_sk_kyber(sk, &skpv);
  pack_pk_kyber(pk, &pkpv, publicseed);
}


/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - unsigned char *c:          pointer to output ciphertext
*              - const unsigned char *m:    pointer to input message (of length KYBER_SYMBYTES bytes)
*              - const unsigned char *pk:   pointer to input public key
*              - const unsigned char *coin: pointer to input random coins used as seed
*                                           to deterministically generate all randomness
**************************************************/
void indcpa_enc(unsigned char *c,
               const unsigned char *m,
               const unsigned char *pk,
               const unsigned char *coins)
{
  polyvec sp, pkpv, ep, at[KYBER_K], bp;
  poly_kyber v, k, epp;
  unsigned char seed[KYBER_SYMBYTES];
  int i;
  unsigned char nonce=0;


  unpack_pk_kyber(&pkpv, seed, pk);

  poly_frommsg(&k, m);

  polyvec_ntt(&pkpv);

  gen_at(at, seed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise(sp.vec+i,coins,nonce++);

  polyvec_ntt(&sp);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise(ep.vec+i,coins,nonce++);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_pointwise_acc(&bp.vec[i],&sp,at+i);

  polyvec_invntt(&bp);
  polyvec_add(&bp, &bp, &ep);

  polyvec_pointwise_acc(&v, &pkpv, &sp);
  poly_invntt_kyber(&v);

  poly_getnoise(&epp,coins,nonce++);

  poly_add_kyber(&v, &v, &epp);
  poly_add_kyber(&v, &v, &k);

  pack_ciphertext(c, &bp, &v);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - unsigned char *m:        pointer to output decrypted message
*              - const unsigned char *c:  pointer to input ciphertext
*              - const unsigned char *sk: pointer to input secret key
**************************************************/
void indcpa_dec(unsigned char *m,
               const unsigned char *c,
               const unsigned char *sk)
{
  polyvec bp, skpv;
  poly_kyber v, mp;

  unpack_ciphertext(&bp, &v, c);
  unpack_sk_kyber(&skpv, sk);

  polyvec_ntt(&bp);

  polyvec_pointwise_acc(&mp,&skpv,&bp);
  poly_invntt_kyber(&mp);

  poly_sub_kyber(&mp, &mp, &v);

  poly_tomsg(m, &mp);
}
