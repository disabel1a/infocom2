#ifndef SETTINGS
#define SETTINGS

/////////////////////////_HAMMING_/////////////////////////

constexpr auto hamming_n = 7;
constexpr auto hamming_k = 4;
constexpr auto hamming_r = hamming_n - hamming_k;
constexpr auto hamming_codes = (1 << hamming_k);

/////////////////////////_CRC_/////////////////////////

constexpr auto crc_r = 16;
constexpr auto crc_m = 8;
constexpr auto crc_k = crc_r + crc_m;
constexpr auto crc_codes = (1 << crc_m);
constexpr auto crc_poly = "10011110101100101";

///////////////////////////////////////////////////////

constexpr auto channels = crc_k / hamming_k;

#endif // !SETTINGS