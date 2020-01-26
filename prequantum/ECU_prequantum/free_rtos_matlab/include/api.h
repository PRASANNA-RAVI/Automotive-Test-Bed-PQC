#ifndef API_H
#define API_H

#include "params.h"

#define CRYPTO_PUBLICKEYBYTES 1472U
#define CRYPTO_SECRETKEYBYTES 3504U
#define CRYPTO_BYTES 2701U

#define CRYPTO_SIGN_ALGNAME "Dilithium_recommended"

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *msg, unsigned long long len,
                const unsigned char *sk);

int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk);


#define CRYPTO_KEM_SECRETKEYBYTES  KYBER_SECRETKEYBYTES
#define CRYPTO_KEM_PUBLICKEYBYTES  KYBER_PUBLICKEYBYTES
#define CRYPTO_KEM_CIPHERTEXTBYTES KYBER_CIPHERTEXTBYTES
#define CRYPTO_KEM_BYTES           KYBER_SYMBYTES

#if   (KYBER_K == 2)
#define CRYPTO_KEM_ALGNAME "Kyber512"
#elif (KYBER_K == 3)
#define CRYPTO_KEM_ALGNAME "Kyber768"
#elif (KYBER_K == 4)
#define CRYPTO_KEM_ALGNAME "Kyber1024"
#else
#error "KYBER_K must be in {2,3,4}"
#endif

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);

int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);


#endif
