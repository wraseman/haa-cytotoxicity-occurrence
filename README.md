# Evaluating regulatory scenarios to limit U.S. nationwide exposure to cytotoxic HAAs
Code and data to support the research described in Peterson et al. 2023:

Eric S. Peterson, William J. Raseman, Benjamin D. Stanford, Gretchen M. Bruce, Heather Klintworth, and David Reckhow (2023) "Evaluating regulatory scenarios to limit U.S. nationwide exposure to cytotoxic HAAs". AWWA Water Science.

This work was funded by [Water Research Foundation Project 5085](https://www.waterrf.org/research/projects/impact-haloacetic-acid-mcl-revision-dbp-exposure-and-health-risk-reduction).

# Contents
- ```python/```: Python code and associated files. 
  - ```data/```: input data for Python scripts. Contains State Dataset, UCMR 4 Dataset, UCMR 4 FOIA Dataset (species-level data), and HAA chemical properties, cytotoxicity, and grouping information.
  - ```output/```: output data and figures created by Python scripts.
  - ```.ipynb files```: Python code in [Jupyter Notebooks](https://jupyter-notebook.readthedocs.io/en/latest/). 
- ```r/```: R code and associated files.
  - ```Data/```: input data for R scripts.
  - ```Output/```: output data and figures created by R scripts. 
  - ```Scripts/```: R scripts to create figures in manuscript. "01", "02", etc. in the file names represent the order script should be run. 
  - ```Supporting Files/```: HAA chemical properties, cytotoxicity, and grouping information.
  - ```.Rproj file```: if using [R Studio](https://en.wikipedia.org/wiki/RStudio), this file automatically sets your working directory to the location of the .Rproj file which will enable the R scripts to run correctly (see [here](https://bookdown.org/ndphillips/YaRrr/projects-in-rstudio.html) for more details). If this file is not used, the user would have to modify each of the file/directory paths to match those on their local machine.

# Abstract
The U.S. Environmental Protection Agency (EPA) is considering a regulatory revision of the Disinfectant and Disinfection Byproduct Rule (DBPR) with a goal of limiting nationwide exposure to DBPs of emerging health concern. The occurrence of four brominated haloacetic acids (HAAs), which are generally more toxic in in vitro assays than the five currently regulated HAAs and are candidates for future regulation, were surveyed in 4,924 public water systems under EPA’s fourth unregulated contaminant monitoring rule (UCMR4). Using UCMR4 data, this study evaluated the nationwide occurrence of nine HAA species and the potential for two regulatory scenarios (the mass sum of all nine HAA species, HAA9, or just the six brominated HAA species, HAA6Br) to control nationwide exposure to the most toxic HAAs. Neither HAA9 nor HAA6Br approaches were effective for identifying water systems that exhibit high HAA exposure, assessed as additive cytotoxicity, because they are more specific to the HAA species that form at high concentrations rather than the species that are most toxic. However, the effectiveness of HAA6Br is highly sensitive to the relative toxicity of one HAA compound, monobromoacetic acid, which has the highest in vitro toxicity among HAAs but also the lowest occurrence and about which little is known regarding in vivo health risks. In contrast to HAA9, systems with high HAA-associated additive toxicity tend to share similar treatment and disinfectant characteristics as systems with high HAA6Br concentrations. Systems with high source water bromide and total organic carbon were far more likely to use chloramines as a disinfectant residual compared to other systems, were no more likely to adopt organic precursor removal technologies (biofiltration, granular activated carbon, and ion exchange) than other systems, on average.

# Partners
This project was completed researchers at Hazen and Sawyer, Intertox, and the University of Massachusetts Amherst.