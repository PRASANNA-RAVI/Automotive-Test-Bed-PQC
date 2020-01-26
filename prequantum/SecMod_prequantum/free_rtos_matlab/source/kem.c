#include "api.h"
// #include "randombytes.h"
#include "fips202.h"
#include "params.h"
#include "verify.h"
#include "indcpa.h"

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
    static unsigned int lfsr_1 = 0x66774433;
    static unsigned int lfsr_2 = 0x0044AABB;
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
* Name:        crypto_kem_keypair
*
* Description: Generates public and private key
*              for CCA-secure Kyber key encapsulation mechanism
*
* Arguments:   - unsigned char *pk: pointer to output public key (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
*              - unsigned char *sk: pointer to output private key (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{
  size_t i;
  indcpa_keypair(pk, sk);
  for(i=0;i<KYBER_INDCPA_PUBLICKEYBYTES;i++)
    sk[i+KYBER_INDCPA_SECRETKEYBYTES] = pk[i];
  sha3_256(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES,pk,KYBER_PUBLICKEYBYTES);
  // for(i=0;i<KYBER_SYMBYTES;i++)
  //   sk[i+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES] = 0xAA;
  // randombytes(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES,KYBER_SYMBYTES);         /* Value z for pseudo-random output on reject */
  get_random_bits(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES,KYBER_SYMBYTES);

  // printf("Public key\n");
  // for(i=0;i<KYBER_INDCPA_PUBLICKEYBYTES;i++)
  // {
  //     printf("%d, ",*(pk+i));
  // }
  // printf("\n");
  //
  // printf("Secret key\n");
  // for(i=0;i<KYBER_SECRETKEYBYTES+2*(KYBER_SYMBYTES);i++)
  // {
  //     printf("%d, ",*(sk+i));
  // }
  // printf("\n");

  return 0;
}

/*************************************************
* Name:        crypto_kem_enc
*
* Description: Generates cipher text and shared
*              secret for given public key
*
* Arguments:   - unsigned char *ct:       pointer to output cipher text (an already allocated array of CRYPTO_CIPHERTEXTBYTES bytes)
*              - unsigned char *ss:       pointer to output shared secret (an already allocated array of CRYPTO_BYTES bytes)
*              - const unsigned char *pk: pointer to input public key (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{
  unsigned char  kr[2*KYBER_SYMBYTES];                                        /* Will contain key, coins */
  unsigned char buf[2*KYBER_SYMBYTES];

  int i;
  // randombytes(buf, KYBER_SYMBYTES);
  get_random_bits(buf, KYBER_SYMBYTES);
  // for(i=0;i<KYBER_SYMBYTES;i++)
  //   buf[i] = 0xAA;

  sha3_256(buf,buf,KYBER_SYMBYTES);                                           /* Don't release system RNG output */

  sha3_256(buf+KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES);                     /* Multitarget countermeasure for coins + contributory KEM */
  sha3_512(kr, buf, 2*KYBER_SYMBYTES);

  indcpa_enc(ct, buf, pk, kr+KYBER_SYMBYTES);                                 /* coins are in kr+KYBER_SYMBYTES */

  sha3_256(kr+KYBER_SYMBYTES, ct, KYBER_CIPHERTEXTBYTES);                     /* overwrite coins in kr with H(c) */
  sha3_256(ss, kr, 2*KYBER_SYMBYTES);                                         /* hash concatenation of pre-k and H(c) to k */
  return 0;
}

/*************************************************
* Name:        crypto_kem_dec
*
* Description: Generates shared secret for given
*              cipher text and private key
*
* Arguments:   - unsigned char *ss:       pointer to output shared secret (an already allocated array of CRYPTO_BYTES bytes)
*              - const unsigned char *ct: pointer to input cipher text (an already allocated array of CRYPTO_CIPHERTEXTBYTES bytes)
*              - const unsigned char *sk: pointer to input private key (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 for sucess or -1 for failure
*
* On failure, ss will contain a randomized value.
**************************************************/
int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{
  size_t i;
  int fail;
  unsigned char cmp[KYBER_CIPHERTEXTBYTES];
  unsigned char buf[2*KYBER_SYMBYTES];
  unsigned char kr[2*KYBER_SYMBYTES];                                         /* Will contain key, coins, qrom-hash */
  const unsigned char *pk = sk+KYBER_INDCPA_SECRETKEYBYTES;

  indcpa_dec(buf, ct, sk);

  for(i=0;i<KYBER_SYMBYTES;i++)                                               /* Multitarget countermeasure for coins + contributory KEM */
    buf[KYBER_SYMBYTES+i] = sk[KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES+i];      /* Save hash by storing H(pk) in sk */
  sha3_512(kr, buf, 2*KYBER_SYMBYTES);

  indcpa_enc(cmp, buf, pk, kr+KYBER_SYMBYTES);                                /* coins are in kr+KYBER_SYMBYTES */

  fail = verify(ct, cmp, KYBER_CIPHERTEXTBYTES);

  sha3_256(kr+KYBER_SYMBYTES, ct, KYBER_CIPHERTEXTBYTES);                     /* overwrite coins in kr with H(c)  */

  cmov(kr, sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES, fail);     /* Overwrite pre-k with z on re-encryption failure */

  sha3_256(ss, kr, 2*KYBER_SYMBYTES);                                         /* hash concatenation of pre-k and H(c) to k */

  return -fail;
}
