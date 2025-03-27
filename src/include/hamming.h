#ifndef HAMMING
#define HAMMING

#include "settings.h"

#include <bitset>

class hamming
{
private:
    uint64_t d;
    std::bitset<hamming_n> gen_matrix[hamming_k];
    std::bitset<hamming_n> check_matrix[hamming_r];
    std::bitset<hamming_r> t_check_matrix[hamming_n];

    std::bitset<hamming_n> code_book[hamming_codes];

    std::bitset<hamming_k> next_value(const std::bitset<hamming_k> &val);

    uint64_t count_weight(std::bitset<hamming_r> &val);

    void init_gen_matrix();
    void init_t_check_matrix();

public:
    hamming();
    ~hamming();

    std::bitset<hamming_n> code(const std::bitset<hamming_k> &msg);
    std::bitset<hamming_n> code(uint64_t idx);

    std::bitset<hamming_k> decode(std::bitset<hamming_n> code);

    void print_hamming_info();
};

#endif // !HAMMING