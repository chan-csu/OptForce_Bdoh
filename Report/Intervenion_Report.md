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


| Reference	| qCO (mmol/gCDW/h) | qCO2 (mmol/gCDW/h) | qH2 (mmol/gCDW/h) |  qAcet (mmol/gCDW/h)	| qetoh (mmol/gCDW/h) | qBDO (mmol/gCDW/h)| u(h-1) | 
| :-------: | :-------: | :-------: | :-------: | :-------: | :-------: | :-------: | :-----:|
| Valgepea 2018 |	-20.042 |	2.125 |	-33.042 |	1.083 |	9.042	| 0 |	0.04 |



Metaclau:


```
load('EXP_Struct.mat')

for i=1:size(EXP_Struct,2)
    Uptakes.rxns{i}=EXP_Struct(i).Name;
    Uptakes.values{i}=EXP_Struct(i).EXP1;                                     %If you want another experiment, add that to EXP_Struct, and change EXP1 accordingly in this line
    Uptakes.type{i}='b';
    metaclau=changeRxnBounds(metaclau,Uptakes.rxns{i},Uptakes.values{i},Uptakes.type{i})
end



Indices=[789,788,790,771,786,787,849];

Biomass_ID=849;

Solution=optimizeCbModel(metaclau,'max','one');
WT_exports=Solution.x(Indices)

>> table(metaclau.rxns(Indices),WT_exports)

ans =

  7×2 table

              Var1              WT_exports
    ________________________    __________

    {'EX_CARBON-MONOXIDE'  }       -20.042
    {'EX_CARBON-DIOXIDE'   }    1.2115e-12
    {'EX_HYDROGEN-MOLECULE'}       -33.042
    {'EX_ACET'             }        3.0048
    {'EX_ETOH'             }        6.0176
    {'EX_BUTANEDIOL'       }             0
    {'EX_BIOMASS'          }      0.086322
    
```

Which is in good agreement with experimental data!


Kinetic model

```
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [1];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [1];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);
        
Kin_Model_Indices=[1 3 2 58 54 62 64 ];
T=array2table(solutions(Kin_Model_Indices,:));
T.Properties.RowNames=model.rxns(Kin_Model_Indices);
    
ans =

  7×18 table

                    Var1        Var2        Var3        Var4        Var5        Var6        Var7        Var8        Var9       Var10       Var11       Var12       Var13       Var14       Var15       Var16       Var17       Var18  
                  ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________

    EX_co          -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042
    EX_co2           4.188        4.86      1.9064      4.1664       3.201      4.4866      4.8634      8.3754      3.3362      3.5127      3.2104      3.7153      3.6796      1.8148      4.8287      5.8895      4.2922      2.8468
    EX_h2          -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042
    EX_ac           5.4151       3.732      5.1397      4.0511      3.8462      5.3852      3.6846      2.6837      4.5681      3.3987      5.5992      4.4483      4.6131      5.0716      4.3638      3.3491      3.7619        4.58
    EX_etoh         1.8208      2.7887      2.9191      2.4959      3.5541      1.4314      2.7828      3.2735      2.2143      3.6919      1.8893      2.6541      2.3637      2.8485      2.2965      2.6636       3.069      3.0595
    EX_bdoh       0.022271    0.093783     0.01219    0.079041     0.01927    0.022496    0.057322     0.04116      0.0198     0.11359    0.050271     0.03696    0.057796    0.085785    0.087939    0.031734    0.059497    0.069689
    EX_biomass    0.014627    0.026152    0.024887    0.024422    0.018093     0.01412    0.017076    0.024029    0.029121    0.016412    0.021893    0.025493    0.024326    0.028383     0.01999    0.027602    0.019099     0.02351               
        
        
        
```
Both models provide reasonably well approximations. 

The following block of code gives the highest bdoh production predicted by FBA:

```
MT=changeObjective(metaclau,'EX_BUTANEDIOL')
solution=optimizeCbModel(MT,'max','one');
MT_exports=solution.x(metaclau_Indices)
table(metaclau.rxns(metaclau_Indices),MT_exports)

ans =

  7×2 table

              Var1              MT_exports
    ________________________    __________

    {'EX_CARBON-MONOXIDE'  }     -20.042  
    {'EX_CARBON-DIOXIDE'   }      1.0578  
    {'EX_HYDROGEN-MOLECULE'}     -33.042  
    {'EX_ACET'             }           0  
    {'EX_ETOH'             }           0  
    {'EX_BUTANEDIOL'       }      4.9485  
    {'EX_BIOMASS'          }           0  

```
This is the highest theoretical bdoh production with the same uptake rates. 
The outputs of the pipeline includes First, Second, and Third order interventions results on the kinetic model. This indicates that only some of the results are significantly higher than wild type bdoh production. Also, some of the bdoh results are lost in the process of mapping simply because that reaction does not exist in the kinetic model. 

### All interventions in on sight

Since there are 18 different parametere sets in the kinetic model, ANOVA is provide a more reasonable comparison of the outputs of the kinetic model. The following script generates all the interventions that were mapped to the kinetic model:

```
load('./Results/First_Order_Core.mat')
%selection of second order OptForce rxns that are in core model
load('./Results/Second_Order_Core.mat')
%selection of third order OptForce rxns that are in core model
load('./Results/Third_Order_Core.mat')

load('./Results/First_Order_results.mat')
%selection of second order OptForce rxns that are in core model
load('./Results/Second_Order_results.mat')
%selection of third order OptForce rxns that are in core model
load('./Results/Third_Order_results.mat')
for i=1:3
    
    First_order_results(:,27,i)= Base_BDOH';
    Second_order_results(:,27,i)=Base_BDOH';
    Third_order_results(:,27,i)=Base_BDOH';

end

```
For First order interventions: 

** Heatmaps

```

First_order_Rxn = model.rxns(First_Order_Core(:,1));
First_order_Intervent = First_Order_Core(:,2)';
First_order_Rxn = [First_order_Rxn;{'WT'}]
First_order_Intervent = [First_order_Intervent 0]
%clear 
clf
subplot(2,1,1)
h1 = heatmap(First_order_results(:,:,3));
h1.XDisplayLabels = First_order_Rxn;
h1ip = get(h1,'InnerPosition');

%Plot the type of intevention of each reaction suggested by Optforce
subplot(2,1,2)
h2 = heatmap(First_order_Intervent);
h2.XDisplayLabels = First_order_Rxn;
set(h2, 'InnerPosition', [h1ip(1) 0.45 h1ip(3)-0.055 0.025]); 

```



