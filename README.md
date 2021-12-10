Parsa Ghadermazi  \
Parsa96@colostate.edu \
1/16/2021

*This project started May-2020, but there has been changes to the structure of
the files. This ReadMe file reflects the latest changes.*

# Introduction

This repository contains the work for predicting effective strategies for 2,3-Butanediol
over production in C. autoethanogenum. First, a genome scale metabolic model is used to predict first, second, and third order strategies for 2,3-Butanediol over production. Then, these predictions are assessed using the existing kinetic ensemble model of C. autoethanogenum [1] to determine the most promising enzymes to target for improving 2,3-BDO production. "Final_Script.m" includes the whole pipeline, and it calls different functions, including the kinetic ensemble model. To reproduce all the results run "Final_Script.m".

This will first run optForce to find the interventions. Then maps the interventions to the core network. Finally, it runs these interventions through the kinetic model by [1]. The intermediates and the results are saved in "./Results/"
To justify using the Metaclau [2] metabolic reconstruction, two GEMMs (Norman et al. 2019; Marcellin et al. 2016) were run on a test, and Metaclau did better on almost all the tests. The script and results of this analysis is reported in "./Report/"


To justify using metaclau, two GEMMs were run on a test, and metaclau did better on almost all the tests. The script and results of this analysis is reported in "./Report/"


./Results/ include the following:
  
- First, Second, and Third_Order_Results. This is the result of kinetic model with three level of perturbation.
- First, Second, and Third_Order_Core. Which include information about the mapped intervention containing the index of the target reaction(s) and their type of regulation (Knockout-0, Upregulation-2, Downregulation-1).
- Base_bdoh is the flux in wildtype by kinetic model for all 18 set of parameters
- Latest_Workspace includes the variables in the workspace after running Final_Script so you don't have to wait to reproduce the results.
- Normalized_First, Second, and Third contains the kinetic model outputs for all interventions normalized by wildtype fluxes.
- MustL, MustLL, MustU, MustUL, MustUU are the outputs of optForce.

**In Response to the reviewers two more analysis were added to the paper:**

1- How would the flux distribution change by changing H2 uptake and holding CO uptake constant.(/Results/H2_Var/)

2- How would the flux distribution change if higher perturbation levels are tested.(/Results/Higher_Pert_Results/)

Both analyses are included in the results directory

# References

[1]- Greene Jennifer, Daniell James, Kopke Michael, Broadbelt Linda, Tyo Keith E.J. (2019). Kinetic ensemble model of gas fermenting Clostridium autoethanogenum for improved ethanol production. Biochemical Engineering Journal. 148: 46-56. DOI: 10.1016/j.bej.2019.04.021.

2- Marcellin Esteban, Behrendorff James B., Nagaraju Shilpa, DeTissera Sashini, Segovia Simon et al. (2016). Low carbon fuels and commodity chemicals from waste gases – systematic approach to understand energy metabolism in a model acetogen. Green Chem. 18(10): 3020-3028. DOI: 10.1039/C5GC02708J. 


