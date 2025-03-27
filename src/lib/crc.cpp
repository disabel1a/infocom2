#include "../include/crc.h"

#include <cmath>
#include <iostream>

crc::crc() :poly_string(crc_poly), d(crc_k), distances(crc_k + 1)
{
    polynome = std::bitset<crc_k>(poly_string + std::string(crc_k - poly_string.size(), '0'));

    std::bitset<crc_m> msg;
    std::bitset<crc_k> code_word;

    uint64_t weight(0);

    for (uint64_t i(0); i < crc_codes; ++i)
    {
        msg = std::bitset<crc_m>(i);
        code_word = code(msg);

        weight = code_word.count();
        ++distances[weight];

        if ((i != 0) && (d > weight))
        {
            d = weight;
        }

        code_book[i] = code_word;
    }
}

crc::~crc()
{
}

std::bitset<crc_k> crc::get_check_part(const std::bitset<crc_m> &msg)
{
    std::bitset<crc_k> remainder(msg.to_string());
    remainder <<= crc_r;

    for (uint64_t i = 0; i < crc_m; ++i) {
        if (remainder.test(crc_k - 1)) {
            remainder ^= polynome;
        }
        remainder <<= 1;
    }

    remainder >>= crc_m;

    return remainder;
}

std::bitset<crc_k> crc::get_syndrom(const std::bitset<crc_k> &code_word)
{
    std::bitset<crc_k> syndrome(code_word); 

    for (uint64_t i = 0; i < crc_m; ++i) {
        if (syndrome.test(crc_k - 1)) {
            syndrome ^= polynome;
        }
        syndrome <<= 1;
    }

    //syndrome >>= crc_k - crc_r;

    return syndrome;
}

uint64_t crc::binomial_coeff(uint64_t n, uint64_t k)
{
    if (k > n)
    {
        return 0;
    }

    if (k == 0 || k == n) 
    {
        return 1;
    }

    uint64_t result = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        result = result * (n - (k - i)) / i;
    }
    return result;
}

std::bitset<crc_k> crc::code(const std::bitset<crc_m> &msg) 
{
    std::bitset<crc_k> output(msg.to_string());
    output <<= crc_r;

    std::bitset<crc_k> remainder = get_check_part(msg);

    output |= remainder;

    return output;
}

std::bitset<crc_k> crc::code(const uint64_t index)
{
    if (index >= crc_codes)
    {
        return code_book[0];
    }

    return code_book[index];
}

std::pair<std::bitset<crc_m>, bool> crc::decode(const std::bitset<crc_k> &code_word) 
{

    std::bitset<crc_m> msg;
    for (int i = 0; i < crc_m; ++i) {
        msg[i] = code_word[i + crc_r];
    }

    bool is_valid = (get_syndrom(code_word) == 0);

    return std::make_pair(msg, is_valid);
}

double crc::upper_border_err(double ch_prob)
{
    double sum = 0.0;

    double ez1 = 0.0;
    double ez2 = 0.0;
    double ez3 = 0.0;

    for (uint64_t i(0); i < d; ++i)
    {
        ez1 = static_cast<double>(binomial_coeff(crc_k, i));
        ez2 = ez1 * std::pow(ch_prob, i);
        ez3 = ez2 * std::pow(1.0 - ch_prob, (crc_k - i));

        sum += ez3;
        //sum += (binomial_coeff(crc_k, i) * std::pow(ch_prob, i) * std::pow(1.0 - ch_prob, (crc_k - i)));
    }

    return 1.0 - sum;
}

double crc::real_err(double ch_prob)
{
    double sum = 0.0;

    double ez1 = 0.0;
    double ez2 = 0.0;

    for (uint64_t i(d); i < crc_k + 1; ++i)
    {
        ez1 = static_cast<double>(distances[i]) * std::pow(ch_prob, i);
        ez2 = ez1 * std::pow(1.0 - ch_prob, (crc_k - i));

        sum += ez2;
        //sum += (distances[i] * std::pow(ch_prob, i) * std::pow(1.0 - ch_prob, (_K_ - i)));
    }

    return sum;
}

void crc::print()
{
    std::cout << "CRC CODES" << std::endl;
    for (uint64_t iter(0); iter < crc_codes; ++iter)
    {
        std::cout << "MSG: " << iter << " CODE: " << code_book[iter] << std::endl;
    }
}
