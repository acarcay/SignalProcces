#ifndef FILTERS_H
#define FILTERS_H

#include <vector>

// --- Moving Average ---
std::vector<double> movingAverage(const std::vector<double>& signal, int windowSize);

// --- 1D Kalman (istenirse hala dursun) ---
class KalmanFilter1D {
private:
    double Q, R, P, x;
public:
    KalmanFilter1D(double processNoise, double measurementNoise, double estimatedError, double initialValue);
    double update(double measurement);
    std::vector<double> apply(const std::vector<double>& signal);
};

// --- 2D Harmonic Oscillator Kalman ---
class KalmanFilter2D {
private:
    double dt;        // time step
    double omega;     // 2*pi*f
    // Process/measurement noise (diagonal Q, scalar R)
    double Qpos, Qvel, R;
    // State [x1; x2] = [position(amplitude); velocity]
    double x1, x2;
    // Covariance P (2x2)
    double P11, P12, P21, P22;

    // Build A for SHO (zero damping)
    inline void A(double& a11, double& a12, double& a21, double& a22) const;

public:
    // Initialize with first measurement as x1, zero velocity by default
    KalmanFilter2D(double dt, double omega,
                   double Qpos, double Qvel, double R,
                   double x1_init, double x2_init = 0.0,
                   double P0_pos = 1.0, double P0_vel = 1.0);

    double update(double z);                       // one-step update with measurement
    std::vector<double> apply(const std::vector<double>& z); // filter entire signal
    double state() const { return x1; }           // current position estimate
};

#endif
std::vector<double> butterworthLowPass(const std::vector<double>& signal,
                                       double dt, double cutoff_hz);

