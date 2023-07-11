---
output:
  pdf_document: default
  html_document: default
---
## Reproducibility Archive

This reproducibility archive allows to reproduce all results & figures of the paper "Time to Intervene: A Continuous-Time Approach to Network Analysis and Centrality". Questions should be addressed to o.ryan@uu.nl. Note that the functions defined in `/functions/` have been updated and compiled in an *R* package `ctnet` available to download from https://github.com/ryanoisin/ctnet/ 

### Functions
The `/functions/` subfolder contains various functions to compute path-specific effects, centrality measures, and to help with plotting

1. `TE.R` Function to compute the CT Total Effect 
2. `DE.R` Function to compute the CT Direct Effect
3. `IE.R` Function to compute the CT Indirect Effect 
4. `ct_centrality.R` Function to compute all CT centrality measures
5. `simPress.R` Function to simulate the effect of an arbitrary press intervention 
6. `EI_VAR.R` Function to compute the one-step and two-step expected influence measures for a DT-VAR model, based on implementations in the `qgraph` package and described in Appendix A.
7. `fun_plots.R` Helper functions to create figures which appear in the main text and process ctsem output

### Empirical Data

Before running `analysis_empirical.R` below, the reader should download the supplementary materials of Kossakowski et al (2017) *Data from ‘critical slowing down as a personalized early warning signal for depression'*. This is available from https://osf.io/c6xt4 and should be placed in the subfolder `/empirical_data/`.

### Main Text and Analysis

1. `analysis_stress.R` Reproduces all figures and tables relating to the toy example, in the sections 1, 2 and 3. Sources all functions above.
2. `analysis_empirical.R` Loads `/empirical_data/ESMdata.csv`. Runs the DT and CT empirical analysis using `ctsem`, saves the output, and creates the figures in Section 4 and tables in Appendix F
3. `SessionInfo.txt` contains information about the R and package versions used for the analysis of empirical data.
4. `estimates.RDS` a file containing the parameter estimates generated by `analysis_empirical.R`
5. `/animations/` files relating to animated pdfs showing how the DT network changes as a funtion of the time-interval
6. `/figures/` contains all figures generated by the code



