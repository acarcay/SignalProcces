#include "filters.h"
#include <cstddef>
#include <cmath>

// -------- Moving Average --------
std::vector<double> movingAverage(const std::vector<double>& signal, int windowSize) {
    if (windowSize < 1) windowSize = 1;
    std::vector<double> filtered(signal.size(), 0.0);
    double sum = 0.0;
    for (size_t i = 0; i < signal.size(); ++i) {
        sum += signal[i];
        if (i >= static_cast<size_t>(windowSize)) {
            sum -= signal[i - windowSize];
            filtered[i] = sum / windowSize;
        } else {
            filtered[i] = sum / (i + 1);
        }
    }
    return filtered;
}

// -------- 1D Kalman (opsiyonel) --------
KalmanFilter1D::KalmanFilter1D(double processNoise, double measurementNoise, double estimatedError, double initialValue)
    : Q(processNoise), R(measurementNoise), P(estimatedError), x(initialValue) {}

double KalmanFilter1D::update(double measurement) {
    P += Q;                              // predict cov
    const double K = P / (P + R);        // gain
    x += K * (measurement - x);          // update state
    P *= (1.0 - K);                      // update cov
    return x;
}

std::vector<double> KalmanFilter1D::apply(const std::vector<double>& signal) {
    std::vector<double> out; out.reserve(signal.size());
    for (double z : signal) out.push_back(update(z));
    return out;
}

// -------- 2D Harmonic Oscillator Kalman --------
inline void KalmanFilter2D::A(double& a11, double& a12, double& a21, double& a22) const {
    const double c = std::cos(omega * dt);
    const double s = std::sin(omega * dt);
    a11 = c;
    a12 = (omega == 0.0) ? dt : (s / omega);
    a21 = -omega * s;
    a22 = c;
}

KalmanFilter2D::KalmanFilter2D(double dt_, double omega_,
                               double Qpos_, double Qvel_, double R_,
                               double x1_init, double x2_init,
                               double P0_pos, double P0_vel)
    : dt(dt_), omega(omega_), Qpos(Qpos_), Qvel(Qvel_), R(R_),
      x1(x1_init), x2(x2_init),
      P11(P0_pos), P12(0.0), P21(0.0), P22(P0_vel) {}

double KalmanFilter2D::update(double z) {
    // ---- Predict ----
    double a11, a12, a21, a22; A(a11,a12,a21,a22);
    const double x1p = a11*x1 + a12*x2;
    const double x2p = a21*x1 + a22*x2;

    // P' = A P A^T + Q  (Q is diagonal)
    double AP11 = a11*P11 + a12*P21;
    double AP12 = a11*P12 + a12*P22;
    double AP21 = a21*P11 + a22*P21;
    double AP22 = a21*P12 + a22*P22;

    double Pp11 = AP11*a11 + AP12*a12 + Qpos;
    double Pp12 = AP11*a21 + AP12*a22;
    double Pp21 = AP21*a11 + AP22*a12;
    double Pp22 = AP21*a21 + AP22*a22 + Qvel;

    // ---- Update ----
    // H = [1 0], so S = Pp11 + R, K = [Pp11; Pp21] / S
    const double S = Pp11 + R;
    const double K1 = Pp11 / S;
    const double K2 = Pp21 / S;

    const double y = z - x1p; // innovation (z - H x_p) since H*x_p = x1p

    x1 = x1p + K1 * y;
    x2 = x2p + K2 * y;

    // P = (I - K H) Pp ; with H=[1 0]
    P11 = (1.0 - K1) * Pp11;
    P12 = (1.0 - K1) * Pp12;
    P21 = Pp21 - K2 * Pp11;
    P22 = Pp22 - K2 * Pp12;

    return x1;
}

// -------- Butterworth Low-Pass (2nd order example) --------
std::vector<double> butterworthLowPass(const std::vector<double>& signal,
                                       double dt, double cutoff_hz) {
    std::vector<double> out(signal.size(), 0.0);
    if (signal.empty()) return out;

    // Örnekleme frekansı ve emniyetli cutoff
    const double fs = 1.0 / dt;
    double fc = std::max(1e-6, std::min(cutoff_hz, 0.499 * fs)); // Nyquist koruması

    // Bilinear transform (pre-warp): K = tan(pi*fc/fs)
    const double K   = std::tan(M_PI * fc / fs);
    const double K2  = K * K;
    const double rt2 = std::sqrt(2.0);

    // Normalize edilmiş katsayılar (a0=1 olacak şekilde)
    const double a0 = 1.0 + rt2 * K + K2;
    const double b0 =  K2 / a0;
    const double b1 =  2.0 * b0;
    const double b2 =  b0;
    const double a1 =  2.0 * (K2 - 1.0) / a0;
    const double a2 =  (1.0 - rt2 * K + K2) / a0;

    // Doğrudan form I
    double x1 = 0.0, x2 = 0.0;
    double y1 = 0.0, y2 = 0.0;
    for (size_t i = 0; i < signal.size(); ++i) {
        const double x0 = signal[i];
        const double y0 = b0*x0 + b1*x1 + b2*x2 - a1*y1 - a2*y2;
        out[i] = y0;
        x2 = x1; x1 = x0;
        y2 = y1; y1 = y0;
    }
    return out;
}



std::vector<double> KalmanFilter2D::apply(const std::vector<double>& z) {
    std::vector<double> out; out.reserve(z.size());
    for (double zi : z) out.push_back(update(zi));
    return out;
}
