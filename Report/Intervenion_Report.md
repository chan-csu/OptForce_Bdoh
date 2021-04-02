# Report on intervention strategies
## Parsa Ghadermazi 4/2/21

-----
The overall pipeline is as follows:

* Using Optforce, top 100 theoretical interventions (First, Second, and Third order interventions) were found. A significant amount of these interventions are not intuitive, if not wrong.

* To have a cross validation tool, an available kinetic model was used. Kinetic models are
ideal tools for examining these interventions since they are directly dealing with enzyme expression levels. Luckily, there was already a kinetic model available for C. autoethanogenum.

* We used the following experimental condition for both metaclau, and the kinetic model:

| Reference  | qCO2 (mmol/gCDW/h) | qH2 (mmol/gCDW/h) |qCO (mmol/gCDW/h)|
| :---------: | :-----------------------------: | :----------: | :-----: |
| Valgepea 2018| 2.125 |	-33.042 | -20.042 |

The reason for choosing this was to be consistent with previously published experimental data, and considering low CO2 emission. 

* The following table shows a comparison between experimentdal data, biomass optimizing metaclau prediction, and the kinetic model predictions: Note that the kinetic model has 18 different parameter sets. We wanted to consider all parameter sets so that the final conclusion is independent of a specific kinetic parameter set. 

Experimental Data:


| Reference	| qCO2 (mmol/gCDW/h) | qCO2 (mmol/gCDW/h) | qH2 (mmol/gCDW/h) |  qAcet (mmol/gCDW/h)	| qetoh (mmol/gCDW/h) | qBDO (mmol/gCDW/h)| u(h-1) | 
| :-------: | :-------: | :-------: | :-------: | :-------: | :-------: | :-------: | :-----:|
| Valgepea 2018 |	-20.042 |	2.125 |	-33.042 |	1.083 |	9.042	| 0 |	0.04 |
