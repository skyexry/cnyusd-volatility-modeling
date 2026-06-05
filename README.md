# CNY/USD Exchange Rate Volatility Modeling

End-to-end ARIMA-GARCH analysis of daily CNY/USD spot exchange rate
(FRED: DEXCHUS, n = 1,248 observations, Nov 2020 – Nov 2025).

## Key Findings

- **Unit root confirmed** (d = 1): ARIMA(0,1,0) selected via AICc minimum (−11310.2) — the series follows a random walk
- **Volatility clustering detected** in squared residuals (Ljung-Box p < 0.001 at lags 12–48), motivating a GARCH extension
- **GARCH(1,1) preferred** over ARCH(1–10): AICc = −11636.60; persistence α₁ + β₁ ≈ 0.982
- **~40% tighter forecast intervals** vs. ARIMA-only: [1.9537, 1.9597] vs. [1.9516, 1.9618]
- Standardized residuals remain leptokurtic (Anderson-Darling p < 0.005); empirical failure rate 6.5% vs. nominal 5%

## Methods

- Stationarity: visual inspection, ACF/PACF, first-difference transformation
- Model selection: AICc grid search over ARIMA(p,1,q), p,q ∈ {0,1,2}
- Heteroskedasticity: Ljung-Box on squared residuals → ARCH(1–10) vs. GARCH(1,1) comparison
- Forecast evaluation: one-step-ahead hold-out on observation 1,249

## Data

Source: [FRED DEXCHUS](https://fred.stlouisfed.org/series/DEXCHUS)
— Chinese Yuan Renminbi to U.S. Dollar, daily, not seasonally adjusted.

## Tools

R · `forecast` · `tseries` (garch) · Minitab (ARIMA selection, residual diagnostics)
