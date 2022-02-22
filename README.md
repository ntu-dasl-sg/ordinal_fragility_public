# Ordinal earthquake fragility curves

This repository contains R code for the analyses in "Order Matters: the Benefits of Ordinal Fragility Curves for Damage and Loss Estimation" by Michele Nguyen and David Lallemant.

- Synthetic_data_experiment.Rmd: The R markdown code for the generation and analysis of the synthetic damage data.
- Nepal_timber_analysis.Rmd: The R markdown code for the case study using damage data from the 2015 Nepal earthquake for buildings with the timber superstructure.   
- KL_x.R: The R code for the computation of the Kullback-Leibler divergence values when the noise in the training data (x%) is varied. 

The output PDFs of the RMarkdown files are also provided for reference. The processed Nepal damage data is available upon request from the authors; however, the raw damage data and ShakeMap is available from [Kathmandu Living Labs](http://eq2015.npc.gov.np/#/) and [United States Geological Survey (USGS)](https://earthquake.usgs.gov/earthquakes/eventpage/us20002926/shakemap/intensity) respectively.
