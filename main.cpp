#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <numbers>

std::pair<double, double> kth_dft(double k, const std::vector<double>& samples) {
    using namespace std::complex_literals;
    std::complex<double> xk = 0.0 + 0.0i;  
    double N = samples.size();  
    for (double n = 0; n < N; ++n) {
        xk += samples[n] * exp((-1.0i) * 2.0 * std::numbers::pi * k * n / N);
    }
    return {xk.real(), xk.imag()};
}

// all of this for rect -> polar basically
double compute_amplitude(double real, double imag, double sample_num) {
    return (2.0 / sample_num) * std::sqrt(real * real + imag * imag); 
    // since dft sums up the amplitude for all samples, 2.0/N makes it so that
    // the actual amplitude is there
}

double compute_phase(double real, double imag) {
    return std::atan2(imag, real);
}

double compute_frequency(int k, int N, double sampling_rate) {
    return k * (sampling_rate / N);
}

int main() {
    double A1, f1, A2, f2;
    std::cout << "input A1: ";
    std::cin >> A1;
    std::cout << "input f1: ";
    std::cin >> f1;
    std::cout << "input A2: ";
    std::cin >> A2;
    std::cout << "input f2: ";
    std::cin >> f2;

    int sampling_time;
    double sampling_rate;
    std::cout << "input sampling time in seconds: ";
    std::cin >> sampling_time;
    std::cout << "input sampling rate (higher = more accurate maybe) in hz: ";
    std::cin >> sampling_rate;

    double sample_num = sampling_rate * sampling_time;

    std::vector<double> samples(sample_num);
    std::vector<double> times(sample_num);

    for (int i = 0; i < sample_num; ++i) {
        times[i] = static_cast<double>(i) / sampling_rate;
    }

    double c1 = 2.0 * std::numbers::pi * f1; // precomputed constants for better performance
    double c2 = 2.0 * std::numbers::pi * f2;
    for (int i = 0; i < sample_num; ++i) {
        double t = times[i];
        samples[i] = A1 * sin(c1 * t) + A2 * sin(c2 * t);
    }   

    std::cout << "\ntime and sample values:\n";
    for (int i = 0; i < sample_num; ++i) {
        double t = times[i];
        std::cout << "t=" << t << " | sample=" << samples[i] << "\n";
    }

    std::cout << "\ndft res:\n";
    for (int k = 0; k < sample_num/2; ++k) {
        std::pair<double, double> xk = kth_dft(static_cast<double>(k), samples);
        double amplitude = compute_amplitude(xk.first, xk.second, sample_num);
        double phase = compute_phase(xk.first, xk.second);
        double frequency = compute_frequency(k, sample_num, sampling_rate);

        if (amplitude > 1e-3) { // very arbitrary, moight change if get feednacl or something
            std::cout << "k=" << k
                      << " | f: " << frequency << " hz"
                      << " | A: " << amplitude
                      << " | phi: " << phase << " rad"
                      << " | eqn: " << amplitude << "*cos(2*pi*" << frequency << "+" << phase << ")\n"; 
        }
    }
    std::cout << "\n+- in output means -";

    return 0;
}