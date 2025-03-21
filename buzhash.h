#ifndef BUZHASH_H
#define BUZHASH_H

#include <stdint.h>


// An implementation of a rolling hash, for fast kmer indexing.
// It's basically a rotate left instruction (top bit to lowest bit)
// followed by and XOR of hash[char] where hash is a randomly generated
// table.

// Precomputed random lookup table.
/*
int main(void) {
    srand(0);
    for (int i = 0; i < 256; i++)
	printf("0x%08x,%c", rand(), " \n"[i%8 == 7]);
    return 0;
}
*/

static uint32_t tab[256] = {
0x6b8b4567, 0x327b23c6, 0x643c9869, 0x66334873,
0x74b0dc51, 0x19495cff, 0x2ae8944a, 0x625558ec,
0x238e1f29, 0x46e87ccd, 0x3d1b58ba, 0x507ed7ab,
0x2eb141f2, 0x41b71efb, 0x79e2a9e3, 0x7545e146,
0x515f007c, 0x5bd062c2, 0x12200854, 0x4db127f8,
0x0216231b, 0x1f16e9e8, 0x1190cde7, 0x66ef438d,
0x140e0f76, 0x3352255a, 0x109cf92e, 0x0ded7263,
0x7fdcc233, 0x1befd79f, 0x41a7c4c9, 0x6b68079a,
0x4e6afb66, 0x25e45d32, 0x519b500d, 0x431bd7b7,
0x3f2dba31, 0x7c83e458, 0x257130a3, 0x62bbd95a,
0x436c6125, 0x628c895d, 0x333ab105, 0x721da317,
0x2443a858, 0x2d1d5ae9, 0x6763845e, 0x75a2a8d4,
0x08edbdab, 0x79838cb2, 0x4353d0cd, 0x0b03e0c6,
0x189a769b, 0x54e49eb4, 0x71f32454, 0x2ca88611,
0x0836c40e, 0x02901d82, 0x3a95f874, 0x08138641,
0x1e7ff521, 0x7c3dbd3d, 0x737b8ddc, 0x6ceaf087,
0x22221a70, 0x4516dde9, 0x3006c83e, 0x614fd4a1,
0x419ac241, 0x5577f8e1, 0x440badfc, 0x05072367,
0x3804823e, 0x77465f01, 0x7724c67e, 0x5c482a97,
0x2463b9ea, 0x5e884adc, 0x51ead36b, 0x2d517796,
0x580bd78f, 0x153ea438, 0x3855585c, 0x70a64e2a,
0x6a2342ec, 0x2a487cb0, 0x1d4ed43b, 0x725a06fb,
0x2cd89a32, 0x57e4ccaf, 0x7a6d8d3c, 0x4b588f54,
0x542289ec, 0x6de91b18, 0x38437fdb, 0x7644a45c,
0x32fff902, 0x684a481a, 0x579478fe, 0x749abb43,
0x3dc240fb, 0x1ba026fa, 0x79a1deaa, 0x75c6c33a,
0x12e685fb, 0x70c6a529, 0x520eedd1, 0x374a3fe6,
0x4f4ef005, 0x23f9c13c, 0x649bb77c, 0x275ac794,
0x39386575, 0x1cf10fd8, 0x180115be, 0x235ba861,
0x47398c89, 0x354fe9f9, 0x15b5af5c, 0x741226bb,
0x0d34b6a8, 0x10233c99, 0x3f6ab60f, 0x61574095,
0x7e0c57b1, 0x77ae35eb, 0x579be4f1, 0x310c50b3,
0x5ff87e05, 0x2f305def, 0x25a70bf7, 0x1dbabf00,
0x4ad084e9, 0x1f48eaa1, 0x1381823a, 0x5db70ae5,
0x100f8fca, 0x6590700b, 0x15014acb, 0x5f5e7fd0,
0x098a3148, 0x799d0247, 0x06b94764, 0x42c296bd,
0x168e121f, 0x1eba5d23, 0x661e3f1e, 0x5dc79ea8,
0x540a471c, 0x7bd3ee7b, 0x51d9c564, 0x613efdc5,
0x0bf72b14, 0x11447b73, 0x42963e5a, 0x0a0382c5,
0x08f2b15e, 0x1a32234b, 0x3b0fd379, 0x68eb2f63,
0x4962813b, 0x60b6df70, 0x06a5ee64, 0x14330624,
0x7fffca11, 0x1a27709e, 0x71ea1109, 0x100f59dc,
0x7fb7e0aa, 0x06eb5bd4, 0x6f6dd9ac, 0x094211f2,
0x00885e1b, 0x76272110, 0x4c04a8af, 0x1716703b,
0x14e17e33, 0x3222e7cd, 0x74de0ee3, 0x68ebc550,
0x2df6d648, 0x46b7d447, 0x4a2ac315, 0x39ee015c,
0x57fc4fbb, 0x0cc1016f, 0x43f18422, 0x60ef0119,
0x26f324ba, 0x7f01579b, 0x49da307d, 0x7055a5f5,
0x5fb8370b, 0x50801ee1, 0x0488ac1a, 0x5fb8011c,
0x6aa78f7f, 0x7672bd23, 0x6fc75af8, 0x6a5f7029,
0x7d5e18f8, 0x5f3534a4, 0x73a1821b, 0x7de67713,
0x555c55b5, 0x3fa62aca, 0x14fce74e, 0x6a3dd3e8,
0x71c91298, 0x09daf632, 0x53299938, 0x1fbfe8e0,
0x5092ca79, 0x1d545c4d, 0x59adea3d, 0x288f1a34,
0x2a155dbc, 0x1d9f6e5f, 0x097e1b4e, 0x51088277,
0x1ca0c5fa, 0x53584bcb, 0x415e286c, 0x7c58fd05,
0x23d86aac, 0x45e6d486, 0x5c10fe21, 0x0e7ffa2b,
0x3c5991aa, 0x4bd8591a, 0x78df6a55, 0x39b7aaa2,
0x2b0d8dbe, 0x6c80ec70, 0x379e21b5, 0x0069e373,
0x2c27173b, 0x4c9b0904, 0x6aa7b75c, 0x1df029d3,
0x5675ff36, 0x3dd15094, 0x3db012b3, 0x2708c9af,
0x5b25ace2, 0x175dfcf0, 0x4f97e3e4, 0x053b0a9e,
0x34fd6b4f, 0x5915ff32, 0x56438d15, 0x519e3149,
0x2c6e4afd, 0x17a1b582, 0x4df72e4e, 0x5046b5a9,
};

static inline uint32_t rotl(uint32_t h) {
    return (h<<1) | (h>>31);
}

static inline uint32_t rotl_k(uint32_t h, uint32_t k) {
    return (h<<k) | (h>>(31-(k-1)));
}

uint32_t hash_init(uint8_t *str, int k) {
    uint32_t h = 0;
    for (int i = 0; i < k; i++)
	h = rotl(h) ^ tab[str[i]];

    return h;
}

static inline uint32_t hash_shift(uint32_t h, uint8_t out, uint8_t in, int k) {
    return rotl(h) ^ rotl_k(tab[out],k) ^ tab[in];
}

#endif // BUZHASH_H
