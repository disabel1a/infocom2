#include "../include/hamming.h"

#include <random>
#include <iostream>

hamming::hamming() : d(hamming_n)
{
    std::bitset<hamming_n> pows(std::string(hamming_k, '0') + std::string(1, '1') + std::string((hamming_r) - 1, '0'));

    for (uint64_t iter(0); iter < hamming_r; ++iter)
    {
        check_matrix[iter] ^= pows;
        pows >>= 1;
    }

    std::bitset<hamming_r> value((1 << hamming_r) - 1);

    for (uint64_t iter(hamming_n - 1); iter >= hamming_r; --iter)
    {
        for (uint64_t col(0); col < hamming_r; ++col)
        {
            check_matrix[col][iter] = value[col];
        }

        do
        {
            value = std::bitset<hamming_r>(value.to_ullong() - 1);
        } while (count_weight(value) <= 1);
    }

    init_gen_matrix();
    init_t_check_matrix();

    std::bitset<hamming_n> code_word;
    std::bitset<hamming_k> curr;
    uint64_t weight(0);

    for (uint64_t iter(0); iter < hamming_codes; ++iter)
    {
        code_word.reset();
        curr = std::bitset<hamming_k>(iter);

        for (uint64_t bit(0); bit < hamming_k; ++bit)
        {
            if (curr[bit])
            {
                code_word ^= gen_matrix[bit];
            }
        }

        weight = code_word.count();

        if ((iter != 0) && (d > weight))
        {
            d = weight;
        }

        code_book[iter] = code_word;
    }
}

uint64_t hamming::count_weight(std::bitset<hamming_r> &val)
{
    uint64_t w(0);
    for (uint64_t iter(0); iter < hamming_r; ++iter)
    {
        if(val[iter])
        {
            ++w;
        }
    }
    return w;
}

hamming::~hamming()
{

}

std::bitset<hamming_k> hamming::next_value(const std::bitset<hamming_k> &val)
{
    std::bitset<hamming_k> one(1);
    return (val[hamming_k - 1]) ? (val << 1) ^ one : val << 1;
}

void hamming::init_gen_matrix()
{
    std::bitset<hamming_n> pows(std::string(1, '1') + std::string(hamming_n - 1, '0'));

    for (uint64_t x(0); x < hamming_k; ++x)
    {
        for (uint64_t y(0); y < hamming_r; ++y)
        {
            gen_matrix[x][(hamming_r - 1) - y] = check_matrix[y][(hamming_n - 1) - x];
        }

        gen_matrix[x] ^= pows;
        pows >>= 1;
    }
}

void hamming::init_t_check_matrix()
{
    for (uint64_t x(0); x < hamming_r; ++x)
    {
        for (uint64_t y(0); y < hamming_n; ++y)
        {
            t_check_matrix[(hamming_n - 1) - y][(hamming_r - 1) - x] = check_matrix[x][y];
        }
    }
}

std::bitset<hamming_n> hamming::code(const std::bitset<hamming_k> &msg)
{
    std::bitset<hamming_n> code_word;
    for (uint64_t bit(0); bit < hamming_k; ++bit)
    {
        if (msg[bit])
        {
            code_word ^= gen_matrix[bit];
        }
    } 

    return code_word;
}

std::bitset<hamming_n> hamming::code(uint64_t idx)
{
    return code_book[idx < hamming_codes ? idx : 0];
}

std::bitset<hamming_k> hamming::decode(std::bitset<hamming_n> code)
{
    std::bitset<hamming_r> syndrome;

    for (uint64_t iter(0); iter < hamming_n; ++iter)
    {
        if (code[(hamming_n - 1) - iter])
        {
            syndrome ^= t_check_matrix[iter];
        }
    }

    if (syndrome.any())
    {
        for (uint64_t iter(0); iter < hamming_n; ++iter)
        {
            if (syndrome == t_check_matrix[iter])
            {
                code[(hamming_n - 1) - iter] = !code[(hamming_n - 1) - iter];
                break;
            }
        }
    }

    std::bitset<hamming_k> msg;
    for (uint64_t iter(0); iter < hamming_k; ++iter)
    {
        msg[iter] = code[(hamming_n - 1) - iter];
    }

    return msg;
}

void hamming::print_hamming_info()
{
    std::cout << "Hamming Info" << std::endl << std::endl;

    std::cout << "Gen Matrix:" << std::endl;
    for (uint64_t iter(0); iter < hamming_k; ++iter)
    {
        std::cout << gen_matrix[iter] << std::endl;
    }

    std::cout << std::endl << "Check Matrix:" << std::endl;
    for (uint64_t iter(0); iter < hamming_r; ++iter)
    {
        std::cout << check_matrix[iter] << std::endl;
    }

    std::cout << std::endl << "Trans Check Matrix:" << std::endl;
    for (uint64_t iter(0); iter < hamming_n; ++iter)
    {
        std::cout << t_check_matrix[iter] << std::endl;
    }

    std::cout << std::endl << "Code Book:" << std::endl;
    for (uint64_t iter(0); iter < hamming_codes; ++iter)
    {
        std::cout << code_book[iter] << std::endl;
    }

    std::cout << std::endl << "Min dist: " << d << std::endl;
}