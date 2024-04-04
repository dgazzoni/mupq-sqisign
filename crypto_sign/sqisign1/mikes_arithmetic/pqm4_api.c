// SPDX-License-Identifier: Apache-2.0

#include <api.h>
#include <string.h>

const char kat_pk[CRYPTO_PUBLICKEYBYTES] = {
    0x9C, 0xD1, 0x15, 0x09, 0x55, 0x1D, 0x64, 0x17, 0x07, 0xA4, 0xD8,
    0x96, 0x58, 0x60, 0xCD, 0x0F, 0xD7, 0x82, 0xCC, 0x1C, 0x87, 0x25,
    0xB5, 0x42, 0xC4, 0xDC, 0x78, 0x5D, 0xEE, 0xCA, 0x56, 0x24, 0x2E,
    0xC6, 0x7D, 0x92, 0x45, 0xCD, 0x46, 0x4B, 0x83, 0x85, 0x54, 0xD3,
    0xA7, 0xFA, 0x4F, 0x1D, 0x90, 0xC1, 0x47, 0x36, 0xF8, 0x24, 0x4E,
    0x21, 0x1D, 0x6E, 0x31, 0xBD, 0xB9, 0x8D, 0x50, 0x0C};

const char kat_sk[CRYPTO_SECRETKEYBYTES] = {
    0x9C, 0xD1, 0x15, 0x09, 0x55, 0x1D, 0x64, 0x17, 0x07, 0xA4, 0xD8, 0x96,
    0x58, 0x60, 0xCD, 0x0F, 0xD7, 0x82, 0xCC, 0x1C, 0x87, 0x25, 0xB5, 0x42,
    0xC4, 0xDC, 0x78, 0x5D, 0xEE, 0xCA, 0x56, 0x24, 0x2E, 0xC6, 0x7D, 0x92,
    0x45, 0xCD, 0x46, 0x4B, 0x83, 0x85, 0x54, 0xD3, 0xA7, 0xFA, 0x4F, 0x1D,
    0x90, 0xC1, 0x47, 0x36, 0xF8, 0x24, 0x4E, 0x21, 0x1D, 0x6E, 0x31, 0xBD,
    0xB9, 0x8D, 0x50, 0x0C, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xAF, 0x71,
    0x87, 0xCD, 0xA3, 0x13, 0xF3, 0x7F, 0xF9, 0xF3, 0xFD, 0xAB, 0xC6, 0x25,
    0x48, 0xF9, 0x69, 0xF2, 0x5C, 0xAB, 0x3B, 0xF9, 0xAD, 0x2C, 0x79, 0x07,
    0x22, 0x7C, 0xA6, 0xAD, 0x10, 0xF3, 0x64, 0x9A, 0xD7, 0xF3, 0x68, 0xB6,
    0xCD, 0x97, 0xD7, 0x6F, 0xF5, 0x67, 0x29, 0xE9, 0x03, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xFA, 0xB3, 0x43, 0xC8, 0x0E, 0xF4, 0x16, 0xB8,
    0xBE, 0x87, 0x4E, 0x2F, 0x6D, 0x2D, 0xA7, 0x38, 0x00, 0x81, 0x24, 0x17,
    0x38, 0x44, 0xDB, 0x60, 0x19, 0xB7, 0xB8, 0x23, 0xC4, 0x9E, 0xB7, 0xD1,
    0x1B, 0x46, 0x6C, 0x99, 0xDC, 0xCD, 0xF9, 0xBD, 0x34, 0x1F, 0x73, 0x0E,
    0x92, 0x32, 0xFE, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xEE, 0x8D,
    0x38, 0xC1, 0x1F, 0xF3, 0xB5, 0x0E, 0x8B, 0xED, 0x00, 0x09, 0x48, 0xC0,
    0x60, 0x4D, 0x63, 0xCB, 0xC8, 0x94, 0x3C, 0xAD, 0x5A, 0x0F, 0xCD, 0x1A,
    0xBF, 0x81, 0x30, 0x81, 0xFB, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xBF, 0xC9, 0x2E, 0x58, 0x2C, 0xC3, 0x50, 0x79,
    0xFA, 0x93, 0x64, 0x3B, 0x1F, 0x12, 0xB3, 0xF0, 0x86, 0x43, 0xF6, 0xBB,
    0xAA, 0x1B, 0xB4, 0x3E, 0xBB, 0x60, 0x1D, 0x25, 0x4D, 0x03, 0x0C, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x8E, 0xF4,
    0xD1, 0x0A, 0x6D, 0x33, 0x6B, 0x5F, 0x52, 0xCD, 0x54, 0x9A, 0x2E, 0x38,
    0x48, 0xD1, 0xDA, 0xE6, 0xDF, 0x6D, 0x25, 0x78, 0x1D, 0x0F, 0x22, 0x06,
    0x23, 0x25, 0xD7, 0xB5, 0xCA, 0x17, 0x00, 0xF2, 0xDC, 0xC3, 0x7B, 0x1C,
    0x21, 0x2B, 0x00, 0x3D, 0xD8, 0x58, 0x66, 0x36, 0xBE, 0x3A, 0x87, 0x89,
    0xA4, 0x05, 0x11, 0x67, 0x93, 0xCA, 0xBE, 0x12, 0xD7, 0xE5, 0xB2, 0x6B,
    0xCC, 0x21, 0x7C, 0x1D, 0x15, 0x74, 0x80, 0x36, 0xC1, 0x80, 0x65, 0xAD,
    0x55, 0x16, 0x20, 0xEF, 0xA2, 0x34, 0xF6, 0x54, 0x76, 0xEE, 0x48, 0xC3,
    0xA5, 0x04, 0x81, 0x18, 0x96, 0x63, 0x36, 0x3D, 0x44, 0x11, 0xE7, 0xED,
    0xD0, 0xB9, 0xBA, 0xC2, 0xC9, 0x3B, 0xD5, 0xB3, 0x49, 0x40, 0x96, 0xB7,
    0xD0, 0xEF, 0xF8, 0xF1, 0xF3, 0xE3, 0x48, 0x1F, 0x38, 0x69, 0xBB, 0xBB,
    0x51, 0xE2, 0x4B, 0x16, 0x1C, 0x24, 0x20, 0x12, 0xF7, 0x9E, 0x4A, 0x62,
    0x50, 0xFE, 0x4C, 0xD6, 0x69, 0x9D, 0x7E, 0x38, 0x24, 0xD3, 0x5E, 0x48,
    0x47, 0xE7, 0x43, 0xA0, 0x7B, 0x09, 0x14, 0x04, 0xA6, 0xBD, 0x71, 0xD7,
    0x5D, 0x2D, 0x7D, 0x09, 0x62, 0xC9, 0xF5, 0x1B, 0x1C, 0xFE, 0x30, 0x0A,
    0xA1, 0xFA, 0xD4, 0xC0, 0x13, 0x1E, 0x69, 0x06, 0x6D, 0x83, 0x2E, 0x3D,
    0x5E, 0x46, 0x20, 0x30, 0xDD, 0x90, 0xE7, 0xD2, 0x5A, 0x28, 0x19, 0x07,
    0x47, 0x5E, 0x62, 0xA3, 0xDE, 0x15, 0x71, 0xCA, 0x21, 0x57, 0xC7, 0x1E,
    0xE3, 0x7F, 0xC1, 0x94, 0x5A, 0x84, 0x04, 0x08, 0xE6, 0xC6, 0xDB, 0x4C,
    0x4D, 0x50, 0x6A, 0x27, 0xEB, 0x24, 0xEF, 0x8E, 0x23, 0x96, 0x66, 0xA6,
    0x63, 0xF6, 0x33, 0xA6, 0xF7, 0x66, 0xDF, 0xE1, 0x27, 0xB3, 0x27, 0x80,
    0xAF, 0xC0, 0x6E, 0x2E, 0x5A, 0x86, 0xD8, 0x6F, 0xB9, 0xA2, 0x53, 0xDD,
    0xD5, 0x31, 0xA3, 0x8E, 0x64, 0x16, 0xE4, 0xFD, 0x5F, 0x12, 0x39, 0x99,
    0x57, 0xCC, 0x87, 0x29, 0x24, 0x6C, 0xB1, 0xD5, 0xD5, 0x72, 0x19, 0xA5,
    0xAA, 0x35, 0xBE, 0x69, 0xA8, 0x21, 0x33, 0x8D, 0x0D, 0x14, 0xE5, 0xAA,
    0x83, 0xE6, 0x3F, 0x84, 0x3C, 0x0E, 0x5C, 0x4F, 0xAB, 0x40, 0xFB, 0xAB,
    0x31, 0x5D, 0x01, 0x9A, 0xDC, 0x0D, 0x29, 0xA2, 0x7B, 0xB5, 0xD9, 0x8E,
    0x39, 0x77, 0xC5, 0x16, 0x31, 0x0A, 0xE2, 0x1D, 0x41, 0xBD, 0x0D, 0xAC,
    0xCC, 0x99, 0x40, 0xA4, 0xCF, 0x2A, 0x36, 0x33, 0xE5, 0x3A, 0xF5, 0x34,
    0xA7, 0xA4, 0x6E, 0xA3, 0x9A, 0x2D, 0x75, 0xAB, 0xBB, 0x41, 0x70, 0xA3,
    0x99, 0x0F, 0xCD, 0xA5, 0xF7, 0x97, 0x91, 0x85, 0x25, 0x2A, 0xAB, 0x38,
    0xA9, 0x27, 0x08, 0x6D, 0xFD, 0x3C, 0x0C, 0x96, 0x04, 0x0F, 0xA6, 0xA8,
    0xC8, 0x77, 0x69, 0xD6, 0xDF, 0xD4, 0xBC, 0x47, 0xBE, 0x0D, 0x38, 0x37,
    0x2C, 0xE5, 0x9B, 0xCD, 0x4C, 0xDD, 0xE3, 0x24, 0xDA, 0x73, 0x17, 0x65,
    0xEB, 0x43, 0x0F, 0x85, 0x84, 0xD7, 0x4D, 0x83, 0x98, 0xC7, 0x3C, 0xB8,
    0x75, 0x4E, 0x91, 0x8C, 0x37, 0x07, 0x61, 0x28, 0xC2, 0xC8, 0x92, 0x62,
    0xA1, 0xA2, 0xBE, 0x43, 0x62, 0xCF, 0xF1, 0x95, 0xD4, 0xA2, 0x76, 0xBE,
    0x74, 0x7E, 0x57, 0x27, 0xFF, 0x58, 0xA7, 0xFE, 0xB0, 0xB4, 0xCF, 0x24,
    0x7F, 0x1B};

typedef struct {
    size_t mlen;
    char msg[59];
    size_t smlen;
    char sm[59 + CRYPTO_BYTES];
} SQISign_KAT_t;

const SQISign_KAT_t kat_lvl1[2] = {
    {
        .mlen = 32,
        .msg = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        .smlen = 32 + CRYPTO_BYTES,
        .sm = {0x71, 0xD4, 0xCF, 0x72, 0xD2, 0x81, 0xEA, 0x05, 0xCE, 0x03,
               0x1F, 0x88, 0xBB, 0x1C, 0xC1, 0x87, 0x97, 0x4E, 0x3D, 0x03,
               0x57, 0xC5, 0x18, 0x8B, 0x4C, 0xA6, 0xEF, 0xEA, 0x90, 0x03,
               0x9D, 0xE0, 0x26, 0x0A, 0x82, 0x9F, 0xF4, 0xC1, 0xF3, 0x07,
               0xDB, 0x59, 0xBE, 0x03, 0x93, 0x6E, 0x8F, 0xF3, 0xB2, 0x01,
               0xB5, 0xB7, 0x1C, 0x31, 0xEF, 0x61, 0xA0, 0x6B, 0xD4, 0x05,
               0x8F, 0xCD, 0x19, 0x96, 0xD2, 0x7C, 0x5C, 0x91, 0xB0, 0x07,
               0x5D, 0x61, 0xF4, 0x5C, 0x27, 0x75, 0xDB, 0x93, 0x0B, 0x05,
               0x26, 0x9E, 0xE8, 0xDD, 0xC8, 0x24, 0x6D, 0xB1, 0x20, 0x04,
               0x33, 0xD2, 0x93, 0x3F, 0xA3, 0xA2, 0x94, 0x9A, 0x2A, 0x06,
               0x75, 0x15, 0xB7, 0x99, 0x12, 0xB6, 0x7B, 0xEB, 0xEB, 0x02,
               0x02, 0x39, 0x32, 0x91, 0xEF, 0x91, 0x76, 0xC8, 0xDE, 0x03,
               0x12, 0x26, 0x8A, 0xFC, 0x4D, 0x31, 0x7C, 0x6A, 0xBE, 0x02,
               0x4D, 0x59, 0xFD, 0x9D, 0xF5, 0x4B, 0x1E, 0x72, 0xB2, 0x01,
               0x01, 0x4D, 0x0B, 0xAB, 0x15, 0xCA, 0x71, 0x24, 0x21, 0x4C,
               0x98, 0x9E, 0x7A, 0xF1, 0x26, 0x04, 0x69, 0x08, 0x03, 0x9F,
               0xD0, 0x50, 0x3F, 0x21, 0x03, 0x3F, 0x84, 0x2A, 0x07, 0xFD,
               0x95, 0xE8, 0x9A, 0xAF, 0x0C, 0x12, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
    },
    {
        .mlen = 59,
        .msg = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
        .smlen = 59 + CRYPTO_BYTES,
        .sm = {0xAD, 0xA1, 0x3C, 0xBB, 0x6C, 0x77, 0xC2, 0x00, 0xE9, 0x07, 0xFD,
               0xD5, 0x0C, 0x3D, 0x2F, 0xA7, 0x86, 0x43, 0x0A, 0x04, 0x1A, 0x62,
               0x17, 0xAF, 0x18, 0xCE, 0x42, 0x41, 0x51, 0x03, 0x25, 0x49, 0x4B,
               0x2E, 0xCF, 0x86, 0x18, 0x53, 0x50, 0x04, 0x1E, 0xFC, 0xE6, 0x5E,
               0x6B, 0x59, 0x38, 0xE0, 0xAE, 0x04, 0x80, 0x68, 0x2D, 0xA6, 0x99,
               0xC8, 0x82, 0xB3, 0x2A, 0x06, 0x03, 0x00, 0xFF, 0xDF, 0x0A, 0x66,
               0x52, 0x1E, 0x3E, 0x07, 0xD3, 0x2C, 0x92, 0xD7, 0x51, 0x51, 0xE2,
               0x71, 0xF6, 0x00, 0x59, 0xE0, 0xBA, 0xCF, 0x9A, 0x24, 0x6F, 0xDE,
               0x70, 0x04, 0x92, 0xC1, 0x80, 0xBC, 0xC9, 0x64, 0x10, 0x86, 0x1D,
               0x02, 0xB9, 0x98, 0xBB, 0xD1, 0xAC, 0x82, 0x72, 0xD6, 0xA6, 0x04,
               0x00, 0xDE, 0x5D, 0x5B, 0x97, 0x4E, 0x9E, 0x1E, 0x37, 0x04, 0x45,
               0x1A, 0x1F, 0x83, 0x37, 0x19, 0xB9, 0x39, 0xA9, 0x03, 0xF9, 0x80,
               0x4E, 0xED, 0xA6, 0x0D, 0xD9, 0x2C, 0x30, 0x06, 0x01, 0x37, 0x65,
               0xFC, 0x85, 0x74, 0x23, 0x81, 0x51, 0x87, 0x8D, 0x77, 0x16, 0xAD,
               0x8E, 0x57, 0x44, 0x03, 0x01, 0x24, 0x2C, 0xCC, 0x33, 0x52, 0x4A,
               0x87, 0xEE, 0xB3, 0x01, 0x7B, 0xBF, 0x4B, 0x79, 0x94, 0x3F, 0xA0,
               0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
               0x00, 0x00, 0x00, 0x00, 0x00},
    },
};

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    memcpy(pk, kat_pk, CRYPTO_PUBLICKEYBYTES);
    memcpy(sk, kat_sk, CRYPTO_SECRETKEYBYTES);
}

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk) {
    for (int i = 0; i < sizeof(kat_lvl1) / sizeof(kat_lvl1[0]); i++) {
        if (mlen == kat_lvl1[i].mlen &&
            memcmp(m, kat_lvl1[i].msg, mlen) == 0) {
            memcpy(sm, kat_lvl1[i].sm, kat_lvl1[i].smlen);
            *smlen = kat_lvl1[i].smlen;
            return 0;
        }
    }

    return 1;
}
