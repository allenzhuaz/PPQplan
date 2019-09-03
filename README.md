[![Travis build
status](https://travis-ci.org/allenzhuaz/PPQplan.svg?branch=master)](https://travis-ci.org/allenzhuaz/PPQplan)
[![](https://www.r-pkg.org/badges/version/PPQplan?color=orange)](https://cran.r-project.org/package=PPQplan)
[![](http://cranlogs.r-pkg.org/badges/grand-total/PPQplan?color=blue)](https://cran.r-project.org/package=PPQplan)
[![](https://img.shields.io/badge/lifecycle-stable-freshgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

Overview
========

This package provides several assessment functions for
statistically-based PPQ sampling plan, including calculating the passing
probability, optimizing the baseline and high performance cutoff points,
visualizing the PPQ plan and power dynamically. The analytical idea is
based on the simulation methods from the textbook “Burdick, R. K.,
LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., &
Zhang, L. (2017). Statistical Methods for CMC Applications. In
Statistical Applications for Chemistry, Manufacturing and Controls (CMC)
in the Pharmaceutical Industry (pp.227-250). Springer, Cham.”

Installation
------------

Open R console, install the pacakge directly from
[CRAN](https://cran.r-project.org/package=PPQplan):

    install.packages("PPQplan")
    library(PPQplan)

Or install the development version from GitHub, first make sure to
install the `devtools` package:

    # install.packages("devtools")
    devtools::install_github("allenzhuaz/PPQplan")
