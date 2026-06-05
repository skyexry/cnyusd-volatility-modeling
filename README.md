# CNY/USD Exchange Rate Volatility Modeling
## ARIMA–GARCH Analysis · R & Minitab

[📄 View Report](https://github.com/skyexry/cnyusd-volatility-modeling/blob/main/report.pdf)

End-to-end volatility analysis of the daily Chinese Yuan to U.S. Dollar
spot exchange rate, combining ARIMA mean modeling (Minitab) with
GARCH conditional-variance estimation (R).

---

## Data

**Source:** [FRED DEXCHUS](https://fred.stlouisfed.org/series/DEXCHUS)
— Chinese Yuan Renminbi to U.S. Dollar, daily, not seasonally adjusted.

`DEXCHUS.csv` contains the raw downloaded series. The first **1,248
observations** (through 2025-11-27) form the estimation sample; observation
1,249 (2025-11-28, logRATE = 1.95658) is held out for forecast evaluation.

---

## Repository Contents

| File | Description |
|------|-------------|
| `DEXCHUS.csv` | Raw daily CNY/USD data from FRED |
| `project2.mpx` | Minitab project: data, transformations, ARIMA output |
| `RES.TXT` | ARIMA(0,1,0) residuals exported from Minitab (n = 1,248) |
| `htfile.txt` | GARCH(1,1) fitted conditional variances h_t from R |
| `analysis.R` | R code for ARCH/GARCH fitting and forecast intervals |
| `report.pdf` | Full write-up with plots and derivations |

---

## Methods

### Part 1–3 · ARIMA (Minitab)

- Computed `logRATE` and `dlogRATE`; confirmed unit root visually and
  via ACF/PACF — `dlogRATE` is approximately white noise
- Fit ARIMA(p,1,q) for p,q ∈ {0,1,2} with and without constant;
  selected **ARIMA(0,1,0) without constant** (random walk) by AICc
- Ljung-Box on residuals: no serial correlation in levels, but squared
  residuals show significant autocorrelation → GARCH extension needed

### Part 4–7 · GARCH (R)

- Loaded ARIMA residuals from `RES.TXT`
- Fit ARCH(q) for q = 0,…,10 via `garch(x, order = c(0, q))`
- Compared AICc: best ARCH is ARCH(4), but **GARCH(1,1) dominates**
- Produced 95% one-step-ahead forecast interval using fitted h_{n+1}

### Part 8–10 · Evaluation

- Standardized residuals remain leptokurtic (Anderson–Darling p < 0.005)
- Empirical failure rate: 81/1246 ≈ **6.5%** vs. nominal 5%

---

## Key Results

| | ARIMA only | ARIMA–GARCH(1,1) |
|---|---|---|
| Model | ARIMA(0,1,0) | ARIMA(0,1,0) + GARCH(1,1) |
| AICc | −11310.2 | −11636.6 |
| 95% PI (one-step) | [1.95160, 1.96177] | [1.95368, 1.95969] |
| PI width | 0.01017 | 0.00601 |
| Actual (log 1249) | 1.95658 ✓ | 1.95658 ✓ |

GARCH(1,1) coefficients: α₀ = 1.29×10⁻⁷, α₁ = 0.0883, β₁ = 0.8942
(α₁ + β₁ ≈ 0.982 — strong volatility persistence, covariance-stationary)

The ARIMA–GARCH interval is **~40% tighter** than the ARIMA-only interval,
because the GARCH component conditions on the relatively low volatility at
the forecast origin. Both intervals cover the actual value.

---

## Reproduction

### Prerequisites

- Minitab 21+
- R 4.x · packages: `tseries`

### Step 1 — ARIMA (Minitab)

1. Open `project2.mpx`
2. All transformations, model fits, and diagnostic plots are stored in
   the project
3. To re-export residuals: **Stat → Time Series → ARIMA → Storage →
   Residuals** → save to `RES.TXT`

### Step 2 — GARCH (R)

```r
# analysis.R
library(tseries)

# Load ARIMA residuals from Minitab
x <- scan("RES.TXT")
N <- length(x)  # 1248

# ARCH(0) log-likelihood by hand
ll0   <- -0.5 * N * (1 + log(2 * pi * mean(x^2)))
aicc0 <- -2 * ll0 + 2 * (0 + 1) * N / (N - 0 - 2)

# ARCH(q) for q = 1,...,10
results <- data.frame(q = 0:10, logLik = NA, AICc = NA)
results$logLik[1] <- ll0
results$AICc[1]   <- aicc0

for (q in 1:10) {
  fit  <- garch(x, order = c(0, q))
  ll   <- as.numeric(logLik(fit))
  aicc <- -2 * ll + 2 * (q + 1) * N / (N - q - 2)
  results$logLik[q + 1] <- ll
  results$AICc[q + 1]   <- aicc
}

# GARCH(1,1)
fit_g11  <- garch(x, order = c(1, 1))
ll_g11   <- as.numeric(logLik(fit_g11))
aicc_g11 <- -2 * ll_g11 + 2 * 3 * N / (N - 1 - 2)

# Export conditional variances
ht <- fit_g11$fitted.values[, 1]^2
write(ht, "htfile.txt", ncolumns = 1)

# One-step-ahead forecast interval
coef_g11 <- coef(fit_g11)
omega <- coef_g11["a0"]
alpha <- coef_g11["a1"]
beta  <- coef_g11["b1"]

e_last <- tail(x, 1)
h_last <- tail(ht, 1)
h_next <- omega + alpha * e_last^2 + beta * h_last

f_n1      <- 1.95668          # ARIMA point forecast from Minitab
lower_95  <- f_n1 - 1.96 * sqrt(h_next)
upper_95  <- f_n1 + 1.96 * sqrt(h_next)
```

---

## Tools

R · `tseries` · Minitab 21
