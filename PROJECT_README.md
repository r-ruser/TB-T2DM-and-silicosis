# TB–DM–Silicosis integrated project

This workspace implements the supplied protocol with an R-only analysis backend.

## Directory contract

- `data/`: all supplied and downloaded data, source documents, manifests, and checksums.
- `R/`: numbered R scripts only.
- `result/`: numbered audit reports, tables, models, figures, source data, and logs.
- `results/`: legacy empty directory retained unchanged; new outputs are written to `result/`.

## Run order

```powershell
Rscript R/00_bootstrap.R
Rscript R/01_environment_audit.R
Rscript R/02_download_geo.R --mode=manifest
Rscript R/02_download_geo.R --mode=supp-manifest
Rscript R/03_download_who.R
Rscript R/10_gbd_audit_clean.R
Rscript R/11_gbd_proxy_analysis.R
Rscript R/12_gbd_proxy_figures.R
```

To download GEO series matrices after checking the manifest and disk budget:

```powershell
Rscript R/02_download_geo.R --mode=matrix
```

Raw/supplementary GEO downloads are deliberately a separate explicit mode because the
single-cell archives can be very large:

```powershell
Rscript R/02_download_geo.R --mode=supp
```

For faster R-native concurrent downloads, use:

```powershell
Rscript R/02_download_geo.R --mode=supp-parallel
```

The recommended analysis-ready download avoids redundant archives/formats and uses four
independent R workers with byte-size verification:

```powershell
Rscript R/02b_download_selected_geo.R
```

## Scientific boundary

The supplied GBD exports do not contain all protocol-defined primary metrics (TB ASIR,
T2DM ASPR, and silicosis ASPR). They contain age-standardized incidence and mortality rates
for tuberculosis, type 2 diabetes mellitus, and pneumoconiosis. The current GBD outputs are therefore
labelled proxy/sensitivity analyses and must not be presented as the preregistered primary result.
