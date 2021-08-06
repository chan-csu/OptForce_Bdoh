Parsa Ghadermazi  \
Parsa96@colostate.edu \
1/16/2021

*This project started May-2020, but there has been changes to the structure of
the files. This ReadMe file reflects the latest changes.*

# Introduction

This repository contains the work for predicting effective strategies for 2,3-Butanediol
over production in C. autoethanogenum. First a genome scale metabolic model is used
to predict first, second, and third order strategies for Butanediol over production.
Then this predictions are mapped to the existing kinetic model, Greene et al., to
find the verify the results with and independent method. "Final_Script.m" includes
the whole pipeline, and it calls different functions, including the kinetic model.

To reproduce all the results run "Final_Script.m"

This will first run optForce to find the interventions. Then maps the interventions to the core network. 
Finally, it runs these interventions through the kinetic model by Greene et al. The intermediates and the results are saved in "./Results/"

To justify using metaclau, two GEMMs were run on a test, and metaclau did better on almost all the tests. The script and results of this analysis is reported in "./Report/"

./Results/ include the following:

- First, Second, and Third_Order_Results. This is the result of kinetic model with three level of perturbation.
- First, Second, and Third_Order_Core. Which include information about the mapped intervention containing the index of the target reaction(s) and their type of regulation (Knockout-0, Upregulation-2, Downregulation-1).
- Base_bdoh is the flux in wildtype by kinetic model for all 18 set of parameters
- Latest_Workspace includes the variables in the workspace after running Final_Script so you don't have to wait to reproduce the results.
- Normalized_First, Second, and Third contains the kinetic model outputs for all interventions normalized by wildtype fluxes.
- MustL, MustLL, MustU, MustUL, MustUU are the outputs of optForce.


