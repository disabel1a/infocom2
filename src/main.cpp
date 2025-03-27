#include "include/crc.h"

#include "include/crc.h"
#include "include/channel.h"
#include "include/file_tools.h"

#include "include/hamming.h"

#include <iostream>

////////////////////////////////////////////////

const std::string ch_berr_file = "../data/ch_berr.txt";
const std::string ap_berr_file = "../data/ap_berr.txt";

const std::string derr_file = "../data/derr.txt";
const std::string t_berr_file = "../data/t_berr.txt";
const std::string t_derr_file = "../data/t_derr.txt";
const std::string t_upper_derr_file = "../data/t_upper_derr.txt";

const std::string ch_usage_file = "../data/ch_usage.txt";

const std::string t_ap_derr_file = "../data/t_ap_derr.txt";

////////////////////////////////////////////////

const double stddev = 15.0;
const int seed = 10;
const uint64_t iterations = 1000000;
const double min_db = -10.0;
const double max_db = 0.0;
const double duration = 1.0;
const double width = (min_db < 0.0) ? ((max_db < 0.0) ? -min_db + -max_db : -min_db + max_db) : min_db + max_db;

const uint64_t capacity = static_cast<uint64_t>((1.0 / duration) * width);

////////////////////////////////////////////////

uint64_t real_iterations(0);

////////////////////////////////////////////////

inline void check_file_procedure();

inline uint64_t msg_generator();

inline uint64_t count_ch_errors(const std::bitset<hamming_n> &tx, const std::bitset<hamming_n> & rx);

inline uint64_t count_ap_errors(const std::bitset<crc_k> &tx, const std::bitset<crc_k> &rx);

inline double count_ch_berr(const uint64_t &errors);

inline double count_ap_berr(const uint64_t &errors);

inline double count_derr(const uint64_t &decoder_errors);

std::vector<double> bits_to_signal(std::bitset<hamming_n> &bits);

std::bitset<hamming_n> signal_to_bits(std::vector<double> &signals);

void to_channels(const std::bitset<crc_k> &code_word, std::bitset<hamming_k> ch_code_words[channels]);

void from_channels(std::bitset<crc_k> &demod_word, const std::bitset<hamming_k> ch_demod_words[channels]);

int main(/*int argc, char const *argv[]*/)
{
    check_file_procedure();

    ////////////////////////////////////////////////
    srand(seed);

    hamming hamming_inst;
    hamming_inst.print_hamming_info();

    crc crc_inst;
    crc_inst.print();

    channel channel_inst(0.0, stddev);

    ////////////////////////////////////////////////

    bool channel_status(0);

    std::bitset<crc_m> msg(0);
    uint64_t msg_idx(0);

    std::bitset<crc_k> code_word(0);

    std::bitset<hamming_k> ch_code_words [channels];
    std::bitset<hamming_k> ch_demod_words [channels];

    std::bitset<hamming_n> hamming_code_word;
    std::bitset<hamming_n> hamming_demod_word;

    std::bitset<crc_k> demod_word(0);
    std::pair<std::bitset<crc_m>, bool> decoder_output;

    std::vector<double> signal;

    ////////////////////////////////////////////////

    std::vector<double> berr_ch;
    std::vector<double> berr_ap;
    std::vector<double> derr;

    std::vector<double> t_berr;
    std::vector<double> t_upper_derr;
    std::vector<double> t_derr;

    std::vector<double> t_ap_derr;

    berr_ch.reserve(capacity);
    berr_ap.reserve(capacity);
    derr.reserve(capacity);
    t_ap_derr.reserve(capacity);

    t_berr.reserve(capacity);
    t_upper_derr.reserve(capacity);
    t_derr.reserve(capacity);

    uint64_t err_bits(0);
    uint64_t err_decoder(0);

    double ch_prob(0.0);
    double ap_prob(0.0);

    ////////////////////////////////////////////////

    std::vector<double> ch_usage;
    ch_usage.reserve(capacity);

    double m = (double) crc_m;
    double n = (double) (crc_k + channels * hamming_r); 

    uint64_t crc_failed (0);

    ////////////////////////////////////////////////

    std::cout << "Channel modelling" << std::endl;
    std::cout << "[" << std::endl;

    for (double Ydb(min_db); Ydb < max_db; Ydb += duration)
    {
        err_bits = 0;
        err_decoder = 0;
        crc_failed = 0;

        ch_prob = 0.0;
        ap_prob = 0.0;

        channel_status = 1;

        t_berr.push_back
        (
            channel_inst.set_SNR_dB(Ydb)
        );

        for (real_iterations = 0; (real_iterations < UINT64_MAX) && ((real_iterations < iterations) || (err_decoder < 1)); ++real_iterations)
        {
            if (channel_status)
            {
                msg_idx = msg_generator();
                msg = std::bitset<crc_m>(msg_idx);

                code_word = crc_inst.code(msg_idx);

                to_channels(code_word, ch_code_words);
            }

            /////////////////////////////// CHANNEL PARTITION ///////////////////////////////

            for (uint32_t ch(0); ch < channels; ++ch)
            {

                hamming_code_word = hamming_inst.code(ch_code_words[ch].to_ullong());

                signal = bits_to_signal(hamming_code_word);
                signal = channel_inst.add_noise(signal);

                hamming_demod_word = signal_to_bits(signal);

                ch_demod_words[ch] = hamming_inst.decode(hamming_demod_word);

                err_bits = count_ch_errors(hamming_code_word, hamming_demod_word);
                if (err_bits != 0)
                {
                    ch_prob += count_ch_berr(err_bits);
                }
            }

            from_channels(demod_word, ch_demod_words);

            /////////////////////////////// ////////////////// ///////////////////////////////

            decoder_output = crc_inst.decode(demod_word);
            channel_status = decoder_output.second;

            err_bits = count_ap_errors(code_word, demod_word);
            ap_prob += count_ap_berr(err_bits);

            if (channel_status)
            {
                std::bitset<crc_m> check = msg ^ decoder_output.first;
                if (check.to_ullong() != 0)
                {
                    ++err_decoder;
                }
            }
            else
            {
                ++crc_failed;
            }
        }
        

        berr_ch.push_back(ch_prob / (double) real_iterations);
        berr_ap.push_back(ap_prob / (double) real_iterations);
        
        derr.push_back(count_derr(err_decoder));

        double dividend = m * static_cast<double>(real_iterations - crc_failed);
        double divisor = n * static_cast<double>(real_iterations);

        ch_usage.push_back(dividend / divisor);

        std::cout << Ydb << "dB" << std::endl;
    }


    std::cout << "]" << std::endl;

    for (auto &val : t_berr)
    {
        t_upper_derr.push_back(crc_inst.upper_border_err(val));
        t_derr.push_back(crc_inst.real_err(val));
    }

    for (auto &val : berr_ap)
    {
        t_ap_derr.push_back(crc_inst.real_err(val));
    }

    write_vector_file(ch_berr_file, berr_ch);
    write_vector_file(ap_berr_file, berr_ap);

    write_vector_file(derr_file, derr);
    write_vector_file(t_berr_file, t_berr);
    write_vector_file(t_derr_file, t_derr);
    write_vector_file(t_upper_derr_file, t_upper_derr);

    write_vector_file(ch_usage_file, ch_usage);
    write_vector_file(t_ap_derr_file, t_ap_derr);

    return 0;
}

inline void check_file_procedure()
{
    check_for_file(ch_berr_file.c_str());
    check_for_file(ap_berr_file.c_str());
    check_for_file(derr_file.c_str());
    check_for_file(t_berr_file.c_str());
    check_for_file(t_derr_file.c_str());
    check_for_file(t_upper_derr_file.c_str());
    check_for_file(ch_usage_file.c_str());
}

inline uint64_t msg_generator()
{
    return rand() % (1 << crc_m);
}

inline uint64_t count_ch_errors(const std::bitset<hamming_n> &tx, const std::bitset<hamming_n> & rx)
{
    return (tx ^ rx).count();
}

inline uint64_t count_ap_errors(const std::bitset<crc_k> &tx, const std::bitset<crc_k> &rx)
{
    return (tx ^ rx).count();
}

inline double count_ch_berr(const uint64_t &errors)
{
    return static_cast<double>(errors) / static_cast<double>(hamming_n * channels);
}

inline double count_ap_berr(const uint64_t &errors)
{
    return static_cast<double>(errors) / static_cast<double>(crc_k);
}

inline double count_derr(const uint64_t &decoder_errors)
{
    return (double)decoder_errors / (double)real_iterations;
}

std::vector<double> bits_to_signal(std::bitset<hamming_n> &bits)
{
    std::vector<double> signals;
    signals.reserve(hamming_n);

    for (uint64_t i = 0; i < hamming_n; ++i)
    {
        double val = (bits[i] == 0) ? -1.0 : 1.0;
        signals.push_back(val);
    }

    return signals;
}

std::bitset<hamming_n> signal_to_bits(std::vector<double> &signals)
{
    std::bitset<hamming_n> bits;
    for (uint64_t i = 0; i < hamming_n; ++i)
    {
        bits[i] = (signals[i] < 0.0) ? 0 : 1;
    }

    return bits;
}

void to_channels(const std::bitset<crc_k> &code_word, std::bitset<hamming_k> ch_code_words[channels])
{
    for (uint64_t ch(0); ch < channels; ++ch)
    {
        for (uint64_t idx(0); idx < hamming_k; ++idx)
        {
            ch_code_words[ch][idx] = code_word[(ch * hamming_k) + idx];
        }
    }
}

void from_channels(std::bitset<crc_k> &demod_word, const std::bitset<hamming_k> ch_demod_words[channels])
{
    for (uint64_t ch(0); ch < channels; ++ch)
    {
        for (uint64_t idx(0); idx < hamming_k; ++idx)
        {
            demod_word[(ch * hamming_k) + idx] = ch_demod_words[ch][idx];
        }
    }
}
