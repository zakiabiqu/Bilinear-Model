# Application: Niño 3 sea surface temperature

This folder applies the `bilinear` package to the monthly sea surface temperature (SST) of the Niño 3 region
(ERSST v5, NOAA Climate Prediction Center, January 1950 to June 2026).

The SST is split into a Fourier seasonal component with a linear trend and a stochastic component.
The stochastic component is modelled with a subset bilinear model.
Out of 832 candidate models, BL(2,2,1,1) has the smallest rolling-origin RMSE.
The selected model is used to forecast July to December 2026.

| File | Contents |
|---|---|
| `paper_eng.pdf` | Paper in English |
| `paper_idn.pdf` | Paper in Indonesian |
| `analisis_bilinear_nino3.R` | Complete R script (data, seasonal component, identification, model search, diagnostics, forecasts) |
| `ersst5.nino.mth.91-20.ascii` | Data (NOAA CPC Niño indices based on ERSST v5) |
| `kandidat.zip` | Saved results of the 832 candidate models |

## Reproducing the analysis

```r
# install.packages("remotes")
remotes::install_github("zakiabiqu/Bilinear-Model")
install.packages(c("ggplot2", "tseries"))
```

Put the script, the data file and the unzipped `kandidat` folder in the same working directory, then run
`analisis_bilinear_nino3.R`. The model search reads the saved candidates from `kandidat` and only evaluates
the missing ones, so the full search is not repeated. Without the folder the search runs in parallel and takes
about two hours.
