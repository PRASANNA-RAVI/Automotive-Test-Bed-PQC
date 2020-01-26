#ifndef PARAMS_H
#define PARAMS_H

#ifndef KYBER_K
#define KYBER_K 3 /* Change this for different security strengths */
#endif

/* Don't change parameters below this line */

#define KYBER_N 256
#define KYBER_Q 7681

#if   (KYBER_K == 2) /* Kyber512 */
#define KYBER_ETA 5
#elif (KYBER_K == 3) /* Kyber768 */
#define KYBER_ETA 4
#elif (KYBER_K == 4) /*KYBER1024 */
#define KYBER_ETA 3
#else
#error "KYBER_K must be in {2,3,4}"
#endif

#define KYBER_SYMBYTES 32   /* size in bytes of shared key, hashes, and seeds */

#define KYBER_POLYBYTES              416
#define KYBER_POLYCOMPRESSEDBYTES    96
#define KYBER_POLYVECBYTES           (KYBER_K * KYBER_POLYBYTES)
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 352)

#define KYBER_INDCPA_MSGBYTES       KYBER_SYMBYTES
#define KYBER_INDCPA_PUBLICKEYBYTES (KYBER_POLYVECCOMPRESSEDBYTES + KYBER_SYMBYTES)
#define KYBER_INDCPA_SECRETKEYBYTES (KYBER_POLYVECBYTES)
#define KYBER_INDCPA_BYTES          (KYBER_POLYVECCOMPRESSEDBYTES + KYBER_POLYCOMPRESSEDBYTES)

#define KYBER_PUBLICKEYBYTES  (KYBER_INDCPA_PUBLICKEYBYTES)
#define KYBER_SECRETKEYBYTES  (KYBER_INDCPA_SECRETKEYBYTES +  KYBER_INDCPA_PUBLICKEYBYTES + 2*KYBER_SYMBYTES) /* 32 bytes of additional space to save H(pk) */
#define KYBER_CIPHERTEXTBYTES  KYBER_INDCPA_BYTES

#define MODE 2

#define SEEDBYTES 32U
#define CRHBYTES 48U
#define N 256U
#define Q 8380417U
#define QBITS 23U
#define ROOT_OF_UNITY 1753U
#define D 14U
#define GAMMA1 ((Q - 1U)/16U)
#define GAMMA2 (GAMMA1/2U)
#define ALPHA (2U*GAMMA2)

#if MODE == 0
#define K 3U
#define L 2U
#define ETA 7U
#define SETABITS 4U
#define BETA 375U
#define OMEGA 64U

#elif MODE == 1
#define K 4U
#define L 3U
#define ETA 6U
#define SETABITS 4U
#define BETA 325U
#define OMEGA 80U

#elif MODE == 2
#define K 5U
#define L 4U
#define ETA 5U
#define SETABITS 4U
#define BETA 275U
#define OMEGA 96U

#elif MODE == 3
#define K 6U
#define L 5U
#define ETA 3U
#define SETABITS 3U
#define BETA 175U
#define OMEGA 120U

#endif

#define POL_SIZE_PACKED ((N*QBITS)/8)
#define POLT1_SIZE_PACKED ((N*(QBITS - D))/8)
#define POLT0_SIZE_PACKED ((N*D)/8)
#define POLETA_SIZE_PACKED ((N*SETABITS)/8)
#define POLZ_SIZE_PACKED ((N*(QBITS - 3))/8)
#define POLW1_SIZE_PACKED ((N*4)/8)

#define POLVECK_SIZE_PACKED (K*POL_SIZE_PACKED)
#define POLVECL_SIZE_PACKED (L*POL_SIZE_PACKED)
#define PK_SIZE_PACKED (SEEDBYTES + K*POLT1_SIZE_PACKED)
#define SK_SIZE_PACKED (2*SEEDBYTES + (L + K)*POLETA_SIZE_PACKED + CRHBYTES + K*POLT0_SIZE_PACKED)
#define SIG_SIZE_PACKED (L*POLZ_SIZE_PACKED + (OMEGA + K) + (N/8 + 8))

#endif
