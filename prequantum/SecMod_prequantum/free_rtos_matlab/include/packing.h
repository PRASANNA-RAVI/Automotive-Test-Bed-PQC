#ifndef PACKING_H
#define PACKING_H

#include "polyvec.h"
#include "poly.h"

// void pack_pk(unsigned char pk[PK_SIZE_PACKED],
//              const unsigned char rho[SEEDBYTES], const polyveck *t1);
// void pack_sk(unsigned char pk[PK_SIZE_PACKED],
//              const unsigned char rho[SEEDBYTES],
//              const unsigned char key[SEEDBYTES],
//              const unsigned char tr[CRHBYTES],
//              const polyvecl *s1,
//              const polyveck *s2,
//              const polyveck *t0);
// void pack_sig(unsigned char sig[SIG_SIZE_PACKED],
//               const polyvecl *z, const polyveck *h, const poly *c);
//
// void unpack_pk(unsigned char rho[SEEDBYTES], polyveck *t1,
//                const unsigned char sk[SK_SIZE_PACKED]);
// void unpack_sk(unsigned char rho[SEEDBYTES],
//                unsigned char key[SEEDBYTES],
//                unsigned char tr[CRHBYTES],
//                polyvecl *s1,
//                polyveck *s2,
//                polyveck *t0,
//                const unsigned char sk[SK_SIZE_PACKED]);
// void unpack_sig(polyvecl *z, polyveck *h, poly *c,
//                 const unsigned char sig[SIG_SIZE_PACKED]);

void pack_pk_kyber(unsigned char *r, const polyvec *pk, const unsigned char *seed);
void unpack_pk_kyber(polyvec *pk, unsigned char *seed, const unsigned char *packedpk);
void pack_ciphertext(unsigned char *r, const polyvec *b, const poly_kyber *v);
void unpack_ciphertext(polyvec *b, poly_kyber *v, const unsigned char *c);
void pack_sk_kyber(unsigned char *r, const polyvec *sk);
void unpack_sk_kyber(polyvec *sk, const unsigned char *packedsk);

#endif
