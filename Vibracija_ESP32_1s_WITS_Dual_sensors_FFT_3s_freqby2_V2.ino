#include <Wire.h>
#include <arduinoFFT.h>
#include <math.h>

// ============================================================
//                     PROJECT OVERVIEW
// ============================================================
// ESP32 + 2x ADXL345 (I2C) vibration monitor
// - Samples both sensors at 200 Hz into circular buffers
// - Computes RMS (total) and RMSAC (AC RMS) from |a| = sqrt(x^2+y^2+z^2)
// - Frequency estimation (0.1 Hz):
//     * 2..8 Hz  -> ACF (autocorrelation)  [robust at very low freq]
//     * 8..100 Hz -> FFT + HPS (Harmonic Product Spectrum) [robust vs 2x errors]
// - WITS output once per second:
//     &&
//     2711..2716 (sensor1)
//     2811..2816 (sensor2)
//     !!
//
// NOTE: WITS frequency prints as "271615.4" (no space) with 1 decimal.
// If your parser dislikes '.', tell me and we will output Hz_x10 integer instead.
// ============================================================


// ============================================================
//                    ADXL345 DEFINITIONS
// ============================================================

#define ADXL345_ADDR_1   0x53   // Sensor 1 (SDO -> GND)
#define ADXL345_ADDR_2   0x1D   // Sensor 2 (SDO -> VCC)

// ADXL345 register addresses
#define REG_POWER_CTL    0x2D
#define REG_DATA_FORMAT  0x31
#define REG_BW_RATE      0x2C
#define REG_DATAX0       0x32

// ADXL345 sensitivity in mg/LSB for +/-2g, FULL_RES (~3.9 mg/LSB)
const float ADXL_LSB_TO_MG = 3.9f;

// ============================================================
//              SAMPLING / WINDOW / OUTPUT SETTINGS
// ============================================================

const uint16_t MAX_SAMPLES        = 512;      // ~2.56 s window @ 200 Hz
const double   SAMPLE_FREQ        = 200.0;
const uint32_t SAMPLE_INTERVAL_US = (uint32_t)(1000000.0 / SAMPLE_FREQ);

const uint32_t WITS_PERIOD_MS = 1000;

// RMS gates (mg)
const int FFT_MIN_AC_RMS_MG = 5;    // below -> don't even attempt frequency
const int RMS_GATE_MG       = 8;    // below -> keep last freq (stability)

// Frequency ranges
const double FREQ_MIN_HZ = 2.0;     // allow 2-3 Hz if real
const double FREQ_MAX_HZ = 100.0;

// Hybrid crossover
const double HYBRID_CROSSOVER_HZ = 8.0; // below -> ACF; above -> FFT/HPS

// FFT/HPS settings
const int    HPS_MAX_HARMONIC     = 3;   // 1x*2x*3x
const double HPS_CONFIDENCE_DB    = 6.0; // below -> reject and keep last

// Low-frequency validation (only used for FFT results < 8 Hz; we prefer ACF there anyway)
const double LOW_FREQ_VALIDATE_HZ       = 8.0;
const double LOW_FREQ_EXTRA_CONF_DB     = 8.0; // stricter acceptance under 8 Hz

// Weak-signal jump limiter (helps prevent random jumps)
const int    WEAK_RMSAC_MG        = 12;
const double MAX_JUMP_HZ_IF_WEAK  = 12.0;

// ACF settings (2..8 Hz)
const double ACF_F_MIN_HZ = 2.0;
const double ACF_F_MAX_HZ = 8.0;
// Minimum "normalized correlation" to accept ACF estimate
const double ACF_MIN_NCORR = 0.20;

// ============================================================
//                   FFT BUFFERS AND OBJECT
// ============================================================

double vReal[MAX_SAMPLES];
double vImag[MAX_SAMPLES];

// ArduinoFFT object
ArduinoFFT<double> FFT(vReal, vImag, MAX_SAMPLES, SAMPLE_FREQ);

// HPS array in log-domain (Nyquist bins)
double hpsLog[MAX_SAMPLES / 2];

// ============================================================
//                 CIRCULAR BUFFERS FOR SENSORS
// ============================================================

int16_t x1Buf[MAX_SAMPLES], y1Buf[MAX_SAMPLES], z1Buf[MAX_SAMPLES];
int16_t x2Buf[MAX_SAMPLES], y2Buf[MAX_SAMPLES], z2Buf[MAX_SAMPLES];

int idx1 = 0, count1 = 0;
int idx2 = 0, count2 = 0;

bool has1 = false, has2 = false;

double lastFreq1 = 0.0;
double lastFreq2 = 0.0;

// ============================================================
//                    ADXL345 I2C HELPERS
// ============================================================

void writeRegister(uint8_t addr, uint8_t reg, uint8_t value) {
  Wire.beginTransmission(addr);
  Wire.write(reg);
  Wire.write(value);
  Wire.endTransmission();
}

bool initADXL345(uint8_t addr) {
  // Output data rate 200 Hz
  writeRegister(addr, REG_BW_RATE, 0x0B);

  // FULL_RES, +/-2g
  writeRegister(addr, REG_DATA_FORMAT, 0x08);

  // Measurement mode
  writeRegister(addr, REG_POWER_CTL, 0x08);

  delay(10);

  // Quick ping (read POWER_CTL)
  Wire.beginTransmission(addr);
  Wire.write(REG_POWER_CTL);
  if (Wire.endTransmission(false) != 0) return false;

  Wire.requestFrom((int)addr, 1);
  if (Wire.available() < 1) return false;
  (void)Wire.read();
  return true;
}

bool readAccelRaw(uint8_t addr, int16_t &x, int16_t &y, int16_t &z) {
  Wire.beginTransmission(addr);
  Wire.write(REG_DATAX0);
  if (Wire.endTransmission(false) != 0) return false;

  Wire.requestFrom((int)addr, 6);

  uint32_t start = millis();
  while (Wire.available() < 6) {
    if (millis() - start > 10) return false;
  }

  uint8_t x0 = Wire.read();
  uint8_t x1 = Wire.read();
  uint8_t y0 = Wire.read();
  uint8_t y1 = Wire.read();
  uint8_t z0 = Wire.read();
  uint8_t z1 = Wire.read();

  x = (int16_t)((x1 << 8) | x0);
  y = (int16_t)((y1 << 8) | y0);
  z = (int16_t)((z1 << 8) | z0);

  return true;
}

// ============================================================
//                   RMS COMPUTATION (TOTAL + AC)
// ============================================================
// RMS(total) is computed from |a| = sqrt(x^2+y^2+z^2)
// RMSAC is computed from |a|-mean(|a|) to remove gravity and DC drift.

void computeRMS(int16_t *xBuf, int16_t *yBuf, int16_t *zBuf,
                int idx, int count,
                int &rms, int &rmsAC)
{
  if (count <= 0) { rms = 0; rmsAC = 0; return; }

  int n = (count > (int)MAX_SAMPLES) ? (int)MAX_SAMPLES : count;

  // Mean of magnitude
  double sum = 0.0;
  for (int k = 0; k < n; k++) {
    int i = idx - 1 - k;
    if (i < 0) i += MAX_SAMPLES;

    double x = xBuf[i];
    double y = yBuf[i];
    double z = zBuf[i];
    sum += sqrt(x*x + y*y + z*z);
  }
  double mean = sum / (double)n;

  // RMS and RMSAC
  double sumSq = 0.0;
  double sumSqAC = 0.0;

  for (int k = 0; k < n; k++) {
    int i = idx - 1 - k;
    if (i < 0) i += MAX_SAMPLES;

    double x = xBuf[i];
    double y = yBuf[i];
    double z = zBuf[i];

    double mag = sqrt(x*x + y*y + z*z);
    double v   = mag;
    double vAC = mag - mean;

    sumSq   += v   * v;
    sumSqAC += vAC * vAC;
  }

  rms   = (int)round(sqrt(sumSq   / (double)n));
  rmsAC = (int)round(sqrt(sumSqAC / (double)n));
}

// ============================================================
//                  PEAK / BIN HELPERS
// ============================================================

static inline double parabolicInterp(const double *arr, int i) {
  double a = arr[i - 1];
  double b = arr[i];
  double c = arr[i + 1];
  double denom = (a - 2.0*b + c);
  if (fabs(denom) < 1e-12) return (double)i;
  double delta = 0.5 * (a - c) / denom;
  return (double)i + delta;
}

static inline int hzToBin(double hz) {
  return (int)round(hz * (double)MAX_SAMPLES / SAMPLE_FREQ);
}

static inline double binToHz(double bin) {
  return (bin * SAMPLE_FREQ) / (double)MAX_SAMPLES;
}

// ============================================================
//            FFT + HPS FREQUENCY ESTIMATION (8..100 Hz)
// ============================================================
// Uses magnitude signal and HPS to suppress 2x harmonic errors.
// Returns lastFreq if not confident.

double computeFreqFFT_HPS_Validated(int16_t *xBuf, int16_t *yBuf, int16_t *zBuf,
                                    int idx, int count, double lastFreq,
                                    int rmsAC)
{
  if (count < 64) return lastFreq;

  int n = (count > (int)MAX_SAMPLES) ? (int)MAX_SAMPLES : count;

  // Build magnitude signal
  double sum = 0.0;
  for (int k = 0; k < n; k++) {
    int i = idx - 1 - k;
    if (i < 0) i += MAX_SAMPLES;

    double x = xBuf[i];
    double y = yBuf[i];
    double z = zBuf[i];

    double mag = sqrt(x*x + y*y + z*z);
    vReal[k] = mag;
    vImag[k] = 0.0;
    sum += mag;
  }

  // Remove DC
  double mean = sum / (double)n;
  for (int k = 0; k < n; k++) vReal[k] -= mean;

  // Zero-fill
  for (int k = n; k < (int)MAX_SAMPLES; k++) { vReal[k] = 0.0; vImag[k] = 0.0; }

  // FFT
  FFT.windowing(FFT_WIN_TYP_HAMMING, FFT_FORWARD);
  FFT.compute(FFT_FORWARD);
  FFT.complexToMagnitude(); // vReal -> spectrum magnitude

  // Build HPS in log-domain
  const double eps = 1e-12;
  const int nyq = (int)(MAX_SAMPLES / 2);

  for (int i = 0; i < nyq; i++) hpsLog[i] = -1e30;

  int binMin = hzToBin(FREQ_MIN_HZ);
  int binMax = hzToBin(FREQ_MAX_HZ);

  if (binMin < 1) binMin = 1;
  if (binMax > nyq - 2) binMax = nyq - 2;

  // Fundamental bin must allow harmonics up to HPS_MAX_HARMONIC
  int binMaxFund = (nyq - 1) / HPS_MAX_HARMONIC;
  if (binMax > binMaxFund) binMax = binMaxFund;

  double sumLog = 0.0;
  int cntLog = 0;

  for (int i = binMin; i <= binMax; i++) {
    double acc = 0.0;
    for (int h = 1; h <= HPS_MAX_HARMONIC; h++) {
      int j = i * h;
      acc += log(vReal[j] + eps);
    }
    hpsLog[i] = acc;
    sumLog += acc;
    cntLog++;
  }

  if (cntLog <= 0) return lastFreq;
  double avgLog = sumLog / (double)cntLog;

  // Find peak in HPS
  int bestBin = -1;
  double bestVal = -1e30;

  for (int i = binMin; i <= binMax; i++) {
    double v = hpsLog[i];
    if (v > bestVal) {
      bestVal = v;
      bestBin = i;
    }
  }

  if (bestBin < 2 || bestBin >= binMax) return lastFreq;

  // Confidence in dB-like units
  double confDb = (bestVal - avgLog) * 8.685889638; // 20*log10(e)

  // Strict for low-frequency FFT results (<8 Hz), but normally we will use ACF there anyway
  // Still included as extra safety.
  double refinedBin = (double)bestBin;
  if (bestBin > binMin && bestBin < binMax) refinedBin = parabolicInterp(hpsLog, bestBin);
  double freqHz = binToHz(refinedBin);

  if (freqHz < LOW_FREQ_VALIDATE_HZ) {
    if (confDb < (HPS_CONFIDENCE_DB + LOW_FREQ_EXTRA_CONF_DB)) return lastFreq;
  } else {
    if (confDb < HPS_CONFIDENCE_DB) return lastFreq;
  }

  // Weak-signal jump limiter
  if (rmsAC <= WEAK_RMSAC_MG && lastFreq > 0.0) {
    if (fabs(freqHz - lastFreq) > MAX_JUMP_HZ_IF_WEAK) return lastFreq;
  }

  // Clamp and round to 0.1 Hz
  if (freqHz < FREQ_MIN_HZ) freqHz = FREQ_MIN_HZ;
  if (freqHz > FREQ_MAX_HZ) freqHz = FREQ_MAX_HZ;
  freqHz = round(freqHz * 10.0) / 10.0;

  return freqHz;
}

// ============================================================
//               ACF FREQUENCY ESTIMATION (2..8 Hz)
// ============================================================
// Uses magnitude signal, DC removed, then normalized autocorrelation.
// Returns 0.0 if not confident.

double computeFreqACF_Mag(int16_t *xBuf, int16_t *yBuf, int16_t *zBuf,
                          int idx, int count)
{
  if (count < (int)MAX_SAMPLES) {
    // Prefer full window for low frequencies
    return 0.0;
  }

  // Convert desired frequency band to lag range
  // lag = Fs / f
  int lagMin = (int)round(SAMPLE_FREQ / ACF_F_MAX_HZ); // for 8 Hz -> 25 samples
  int lagMax = (int)round(SAMPLE_FREQ / ACF_F_MIN_HZ); // for 2 Hz -> 100 samples

  if (lagMin < 2) lagMin = 2;
  if (lagMax > (int)MAX_SAMPLES - 2) lagMax = (int)MAX_SAMPLES - 2;

  // Build magnitude array in chronological order (oldest->newest) into vReal
  // and remove mean (DC).
  double sum = 0.0;
  for (int k = 0; k < (int)MAX_SAMPLES; k++) {
    int i = (idx - (int)MAX_SAMPLES + k);
    while (i < 0) i += (int)MAX_SAMPLES;
    i %= (int)MAX_SAMPLES;

    double x = xBuf[i];
    double y = yBuf[i];
    double z = zBuf[i];
    double mag = sqrt(x*x + y*y + z*z);

    vReal[k] = mag;
    sum += mag;
  }

  double mean = sum / (double)MAX_SAMPLES;
  for (int k = 0; k < (int)MAX_SAMPLES; k++) vReal[k] -= mean;

  // Compute signal energy for normalization
  double energy = 0.0;
  for (int k = 0; k < (int)MAX_SAMPLES; k++) energy += vReal[k] * vReal[k];
  if (energy < 1e-9) return 0.0;

  // Normalized autocorrelation
  int bestLag = 0;
  double bestNCorr = -1e30;

  for (int lag = lagMin; lag <= lagMax; lag++) {
    double corr = 0.0;
    for (int k = 0; k < (int)MAX_SAMPLES - lag; k++) {
      corr += vReal[k] * vReal[k + lag];
    }
    // Normalize by energy (approx; good enough for gating)
    double nCorr = corr / energy;

    if (nCorr > bestNCorr) {
      bestNCorr = nCorr;
      bestLag = lag;
    }
  }

  // Accept only if correlation is strong enough
  if (bestLag <= 0 || bestNCorr < ACF_MIN_NCORR) return 0.0;

  double freq = SAMPLE_FREQ / (double)bestLag;

  // Round to 0.1 Hz
  freq = round(freq * 10.0) / 10.0;
  return freq;
}

// ============================================================
//                 HYBRID FREQUENCY ESTIMATION
// ============================================================
// - Gate by RMSAC for stability
// - Prefer FFT/HPS for >=8 Hz
// - Prefer ACF for <8 Hz (only if ACF confident)

double computeFreqHybrid(int16_t *xBuf, int16_t *yBuf, int16_t *zBuf,
                         int idx, int count, double lastFreq,
                         int rmsAC)
{
  // Hard stability gate
  if (rmsAC < RMS_GATE_MG) {
    return lastFreq;
  }

  // Try FFT/HPS first
  double fFFT = computeFreqFFT_HPS_Validated(xBuf, yBuf, zBuf, idx, count, lastFreq, rmsAC);

  if (fFFT >= HYBRID_CROSSOVER_HZ) {
    return fFFT;
  }

  // If FFT suggests low freq, use ACF (more reliable at low Hz)
  double fACF = computeFreqACF_Mag(xBuf, yBuf, zBuf, idx, count);

  if (fACF > 0.0) {
    // Use ACF only inside ACF band; otherwise keep last
    if (fACF >= ACF_F_MIN_HZ && fACF <= ACF_F_MAX_HZ) return fACF;
  }

  // If ACF not confident, keep last (avoids 2-5 Hz noise "stickiness")
  return lastFreq;
}

// ============================================================
//                           SETUP
// ============================================================

void setup() {
  Serial.begin(115200);
  Wire.begin();

  has1 = initADXL345(ADXL345_ADDR_1);
  has2 = initADXL345(ADXL345_ADDR_2);

  delay(100);
}

// ============================================================
//                            LOOP
// ============================================================

void loop() {
  static uint32_t lastSampleUs = 0;
  static uint32_t lastWitsMs   = 0;

  uint32_t nowUs = micros();
  uint32_t nowMs = millis();

  // --------------------------------------------------------
  // 1) Sampling (200 Hz)
  // --------------------------------------------------------
  if (nowUs - lastSampleUs >= SAMPLE_INTERVAL_US) {
    lastSampleUs += SAMPLE_INTERVAL_US;

    int16_t x, y, z;

    // Sensor 1
    if (has1) {
      if (readAccelRaw(ADXL345_ADDR_1, x, y, z)) {
        x1Buf[idx1] = (int16_t)round(x * ADXL_LSB_TO_MG);
        y1Buf[idx1] = (int16_t)round(y * ADXL_LSB_TO_MG);
        z1Buf[idx1] = (int16_t)round(z * ADXL_LSB_TO_MG);

        idx1 = (idx1 + 1) % (int)MAX_SAMPLES;
        if (count1 < (int)MAX_SAMPLES) count1++;
      } else {
        has1 = false;
      }
    }

    // Sensor 2
    if (has2) {
      if (readAccelRaw(ADXL345_ADDR_2, x, y, z)) {
        x2Buf[idx2] = (int16_t)round(x * ADXL_LSB_TO_MG);
        y2Buf[idx2] = (int16_t)round(y * ADXL_LSB_TO_MG);
        z2Buf[idx2] = (int16_t)round(z * ADXL_LSB_TO_MG);

        idx2 = (idx2 + 1) % (int)MAX_SAMPLES;
        if (count2 < (int)MAX_SAMPLES) count2++;
      } else {
        has2 = false;
      }
    }
  }

  // --------------------------------------------------------
  // 2) WITS output (1 Hz)
  // --------------------------------------------------------
  if (nowMs - lastWitsMs >= WITS_PERIOD_MS) {
    lastWitsMs = nowMs;

    Serial.println("&&");

    // =====================================================
    //                SENSOR 1 → 2711–2716
    // =====================================================
    if (has1 && count1 > 0) {
      int i = idx1 - 1;
      if (i < 0) i = MAX_SAMPLES - 1;

      int X = abs(x1Buf[i]);
      int Y = abs(y1Buf[i]);
      int Z = abs(z1Buf[i]);

      int rms, rmsAC;
      computeRMS(x1Buf, y1Buf, z1Buf, idx1, count1, rms, rmsAC);

      double freq = lastFreq1;

      if (rmsAC >= FFT_MIN_AC_RMS_MG) {
        freq = computeFreqHybrid(x1Buf, y1Buf, z1Buf, idx1, count1, lastFreq1, rmsAC);
        lastFreq1 = freq;
      } else {
        // Below minimum: keep last (or set 0.0 if preferred)
        // freq = 0.0;
      }

      Serial.print("2711"); Serial.println(X);
      Serial.print("2712"); Serial.println(Y);
      Serial.print("2713"); Serial.println(Z);
      Serial.print("2714"); Serial.println(rms);
      Serial.print("2715"); Serial.println(rmsAC);
      Serial.print("2716"); Serial.println(freq, 1);
    }

    // =====================================================
    //                SENSOR 2 → 2811–2816
    // =====================================================
    if (has2 && count2 > 0) {
      int i = idx2 - 1;
      if (i < 0) i = MAX_SAMPLES - 1;

      int X = abs(x2Buf[i]);
      int Y = abs(y2Buf[i]);
      int Z = abs(z2Buf[i]);

      int rms, rmsAC;
      computeRMS(x2Buf, y2Buf, z2Buf, idx2, count2, rms, rmsAC);

      double freq = lastFreq2;

      if (rmsAC >= FFT_MIN_AC_RMS_MG) {
        freq = computeFreqHybrid(x2Buf, y2Buf, z2Buf, idx2, count2, lastFreq2, rmsAC);
        lastFreq2 = freq;
      }

      Serial.print("2811"); Serial.println(X);
      Serial.print("2812"); Serial.println(Y);
      Serial.print("2813"); Serial.println(Z);
      Serial.print("2814"); Serial.println(rms);
      Serial.print("2815"); Serial.println(rmsAC);
      Serial.print("2816"); Serial.println(freq, 1);
    }

    Serial.println("!!");
  }
}