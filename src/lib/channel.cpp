#include "../include/channel.h"
#include <cmath>

channel::channel(double mean, double stddev) : mean(mean), stddev(stddev), dist(mean, stddev)
{
}

channel::~channel()
{
}

double channel::set_SNR_dB(double SNR_dB) 
{
    double snr_linear = std::pow(10.0, SNR_dB / 10.0);
    double noise_std = std::sqrt(1.0 / (2.0 * snr_linear));
    reset_params(0.0, noise_std);

    return 0.5 * std::erfc(std::sqrt(snr_linear));
}

void channel::reset_params(double mean, double stddev)
{
    dist = std::normal_distribution<double>(mean, stddev);
}

std::vector<double> channel::add_noise(std::vector<double> &singals)
{
    std::vector<double> output;
    output.reserve(singals.size());
    
    for (auto &s : singals)
    {
        output.push_back(s + dist(generator));
    }

    return output;
}
