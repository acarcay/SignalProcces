# Signal Processing Demo (C++)

This project demonstrates basic and advanced signal processing techniques in C++.

## Features
- **Signal generation**: sinusoidal signal with Gaussian noise
- **Moving Average filter**
- **Kalman filters**:
  - 1D Kalman
  - 2D Harmonic Oscillator Kalman
- **Butterworth low-pass filter (2nd order)**
- **CSV export** for offline analysis
- **Live plotting** using [matplotlib-cpp](https://github.com/lava/matplotlib-cpp)

## Requirements
- C++17
- CMake
- Python 3 with `matplotlib` and `numpy`
- [matplotlib-cpp](https://github.com/lava/matplotlib-cpp)

## Build
```bash
mkdir build && cd build
cmake ..
make
