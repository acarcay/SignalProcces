#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include "filters.h"
#include "matplotlibcpp.h"   // grafik için ekleme

namespace plt = matplotlibcpp;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Params {
    double A = 1.0;          // amplitude
    double f = 1.0;          // frequency (Hz)
    double dt = 0.01;        // step (s)
    double T  = 10.0;        // total time (s)
    double noise_std = 0.2;  // Gaussian noise sigma
    int    windowSize = 31;  // MA window
    // 2D Kalman (harmonic oscillator) params
    double Qpos = 1e-4;
    double Qvel = 1e-2;
    double R    = 0.04;
    double P0pos = 1.0, P0vel = 1.0;
};

static void write_csv(const std::string& path,
                      const std::vector<double>& t,
                      const std::vector<double>& clean,
                      const std::vector<double>& noisy,
                      const std::vector<double>& ma,
                      const std::vector<double>& kal2d)
{
    std::ofstream f(path);
    f << "Time,Clean,Noisy,MA_Filtered,Kalman2D_Filtered\n";
    for (size_t i = 0; i < t.size(); ++i) {
        f << t[i] << "," << clean[i] << "," << noisy[i] << ","
          << ma[i] << "," << kal2d[i] << "\n";
    }
}

int main(int argc, char** argv) {
    Params p;

    // Basit CLI: ./sim [windowSize] [Qpos] [Qvel] [R]
    if (argc >= 2) p.windowSize = std::max(1, std::stoi(argv[1]));
    if (argc >= 3) p.Qpos = std::stod(argv[2]);
    if (argc >= 4) p.Qvel = std::stod(argv[3]);
    if (argc >= 5) p.R    = std::stod(argv[4]);

    // RNG
    std::random_device rd; std::mt19937 gen(rd());
    std::normal_distribution<double> n01(0.0, p.noise_std);

    const double omega = 2.0 * M_PI * p.f;

    // Sinyal
    std::vector<double> time, clean, noisy;
    for (double t = 0.0; t <= p.T + 1e-12; t += p.dt) {
        const double s = p.A * std::sin(omega * t);
        const double z = s + n01(gen);
        time.push_back(t); clean.push_back(s); noisy.push_back(z);
    }

    // Filtreler
    auto ma = movingAverage(noisy, p.windowSize);

    auto butter = butterworthLowPass(noisy, p.dt, 1.5);  // cutoff freq = 2 Hz (örnek)

    KalmanFilter2D kf2(p.dt, omega, p.Qpos, p.Qvel, p.R,
                       noisy.front(), 0.0, p.P0pos, p.P0vel);
    auto kal2d = kf2.apply(noisy);

    // Çıktı CSV
    write_csv("filtered_signals.csv", time, clean, noisy, ma, kal2d);

    {
        std::ofstream f("signal_data.csv");
        f << "Time,Clean,Noisy\n";
        for (size_t i = 0; i < time.size(); ++i)
            f << time[i] << "," << clean[i] << "," << noisy[i] << "\n";
    }

    std::cout << "Saved 'filtered_signals.csv' with 2D Kalman and MA.\n";


    // 📊 Grafik
    plt::figure();
    plt::named_plot("Noisy", time, noisy, "r.");
    plt::named_plot("Clean", time, clean, "b-");
    plt::named_plot("Moving Avg", time, ma, "g-");
    plt::named_plot("Kalman2D", time, kal2d, "k-");
    plt::named_plot("Butterworth", time, butter, "m-");
    plt::legend();
    plt::xlabel("Time (s)");
    plt::ylabel("Amplitude");
    plt::title("Signal Filtering Demo");
    plt::show();


    return 0;
}
