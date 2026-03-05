ESP32 Dual ADXL345 Vibration Monitor

This project measures vibration amplitude and dominant frequency using two ADXL345 accelerometers connected to an ESP32.

The system samples acceleration at 200 Hz, computes RMS vibration level, and estimates the dominant vibration frequency.

Frequency estimation uses a hybrid algorithm to improve stability and avoid harmonic errors.

Low frequencies (2–8 Hz) are detected using autocorrelation (ACF).
Higher frequencies (8–100 Hz) are detected using FFT with Harmonic Product Spectrum (HPS).

This approach provides stable frequency estimation even when vibration amplitude is small.

Main features

- Supports two ADXL345 sensors
- 200 Hz sampling rate
- RMS and AC RMS calculation
- Hybrid frequency detection:
  - ACF for 2–8 Hz
  - FFT + HPS for 8–100 Hz
- Frequency resolution: 0.1 Hz
- Noise filtering using RMS gates
- Stable frequency tracking
- WITS formatted serial output
- Output update every second

Hardware

- ESP32
- 2x ADXL345 accelerometer (I2C)

Sensor I2C addresses

Sensor 1: 0x53  
Sensor 2: 0x1D

Sampling parameters

Sampling frequency: 200 Hz  
FFT window: 512 samples (~2.56 seconds)

Output format (WITS)

Data is transmitted once per second:

&&
2711  X (mg)
2712  Y (mg)
2713  Z (mg)
2714  RMS total (mg)
2715  RMS AC (mg)
2716  Dominant frequency (Hz)

2811  X (mg)
2812  Y (mg)
2813  Z (mg)
2814  RMS total (mg)
2815  RMS AC (mg)
2816  Dominant frequency (Hz)
!!

Notes

Frequency is reported with 0.1 Hz resolution.

RMS AC gating prevents unstable frequency detection when vibration amplitude is very small.

Hybrid ACF + FFT/HPS algorithm prevents common octave errors where FFT detects double frequency instead of the fundamental.
