#ifndef CRC
#define CRC 

#include "settings.h"

#include <bitset>
#include <vector>

class crc
{
private:
    std::string poly_string;
    std::bitset<crc_k> polynome;

    std::bitset<crc_k> code_book[crc_codes];

    uint64_t d;
    std::vector<uint64_t> distances;

    std::bitset<crc_k> get_check_part(const std::bitset<crc_m> &msg);
    std::bitset<crc_k> get_syndrom(const std::bitset<crc_k> &code_word);

    uint64_t binomial_coeff(uint64_t n, uint64_t k);

public:
    crc();
    ~crc();

    std::bitset<crc_k> code(const std::bitset<crc_m> &msg);
    std::bitset<crc_k> code(const uint64_t index);

    std::pair<std::bitset<crc_m>, bool> decode(const std::bitset<crc_k> &code_word);

    double upper_border_err(double ch_prob);
    double real_err(double ch_prob);

    void print();
};

#endif // !CRC