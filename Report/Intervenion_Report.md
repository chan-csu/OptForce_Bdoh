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

**Heatmaps**

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

![First_Order_heatmap](https://github.com/chan-csu/OptForce_Bdoh/blob/master/Report/First_Order_Heatmap.jpg)

From the heatmap 'PFOR' enzyme upregualtion seems like a candidate. To be more certain, the following block provides ANOVA test:

```

[p,tbl,stats] = anova1(First_order_results(:,:,3));
[c,m,h,gnames] = multcompare(stats);

```

![First_Order_ANOVA](https://github.com/chan-csu/OptForce_Bdoh/blob/master/Report/ANOVA_First.jpg)

As it can be seen, PFOR upregulation results in significantly higher bdoh production with blue showing wild type bdoh production.

Same for Second and Third order interventions:

```

Second_order_Rxn_A = model.rxns(Second_Order_Core(:,1));
Second_order_Rxn_A=[Second_order_Rxn_A;{'WT'}]
Second_order_Rxn_B = model.rxns(Second_Order_Core(:,3));
Second_order_Rxn_B=[Second_order_Rxn_B;{'WT'}]
Second_order_Intervent_A = Second_Order_Core(:,2)';
Second_order_Intervent_A=[Second_order_Intervent_A 0];
Second_order_Intervent_B = Second_Order_Core(:,4)';
Second_order_Intervent_B=[Second_order_Intervent_B 0];

subplot(3,1,1)
h1 = heatmap(Second_order_results(:,:,3));
h1ip = get(h1,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(3,1,2)
h2 = heatmap(Second_order_Intervent_A);
h2.ColorLimits = [0 2];
h2.XDisplayLabels = Second_order_Rxn_A;
set(h2, 'InnerPosition', [h1ip(1) 0.65 h1ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(3,1,3)
h3 = heatmap(Second_order_Intervent_B);
h3.ColorLimits = [0 2];
h3.XDisplayLabels = Second_order_Rxn_B;
set(h3, 'InnerPosition', [h1ip(1) 0.55 h1ip(3)-0.055 0.025]);

[p,tbl,stats] = anova1(Second_order_results(:,:,3));
[c,m,h,gnames] = multcompare(stats);

```
![First_Order_ANOVA](https://github.com/chan-csu/OptForce_Bdoh/blob/master/Report/Heatmap_Second.jpg)

![First_Order_ANOVA](https://github.com/chan-csu/OptForce_Bdoh/blob/master/Report/ANOVA_Second.jpg)


Three second order interventions are detected to be significant:

{'ACLDC'} --> Upregulation , {'FBP' } --> Knockout


{'BTDDx'} --> Upregulation,  {'PFOR' } --> Upregulation


{'ACLS' } --> Upregulation,  {'GLUD'} --> Knockout

And, third order interventions:

```

Third_order_Rxn_A = model.rxns(Third_Order_Core(:,1));
Third_order_Rxn_A = [Third_order_Rxn_A;{'WT'}]
Third_order_Rxn_B = model.rxns(Third_Order_Core(:,3));
Third_order_Rxn_B = [Third_order_Rxn_B;{'WT'}]
Third_order_Rxn_C = model.rxns(Third_Order_Core(:,5));
Third_order_Rxn_C = [Third_order_Rxn_C;{'WT'}]
Third_order_Intervent_A = Third_Order_Core(:,2)';
Third_order_Intervent_A =[Third_order_Intervent_A  0]
Third_order_Intervent_B = Third_Order_Core(:,4)';
Third_order_Intervent_B =[Third_order_Intervent_B  0]
Third_order_Intervent_C = Third_Order_Core(:,6)';
Third_order_Intervent_C =[Third_order_Intervent_C  0]

subplot(4,1,1)
h1 = heatmap(Third_order_results(:,:,3));
h1ip = get(h1,'InnerPosition');

%Plot the type of intevention of each reaction A suggested by Optforce
subplot(4,1,2)
h2 = heatmap(Third_order_Intervent_A);
h2.ColorLimits = [0 2];
h2.XDisplayLabels = Third_order_Rxn_A;
set(h2, 'InnerPosition', [h1ip(1) 0.65 h1ip(3)-0.055 0.025]); 

%Plot the type of intevention of each reaction B suggested by Optforce
subplot(4,1,3)
h3 = heatmap(Third_order_Intervent_B);
h3.ColorLimits = [0 2];
h3.XDisplayLabels = Third_order_Rxn_B;
set(h3, 'InnerPosition', [h1ip(1) 0.55 h1ip(3)-0.055 0.025]);

subplot(4,1,4)
h4 = heatmap(Third_order_Intervent_C);
h4.ColorLimits = [0 2];
h4.XDisplayLabels = Third_order_Rxn_C;
set(h4, 'InnerPosition', [h1ip(1) 0.45 h1ip(3)-0.055 0.025]);

[p,tbl,stats] = anova1(Third_order_results(:,:,3));
[c,m,h,gnames] = multcompare(stats);

```
![First_Order_ANOVA](https://github.com/chan-csu/OptForce_Bdoh/blob/master/Report/ANOVA_Third.jpg)

![First_Order_ANOVA](https://github.com/chan-csu/OptForce_Bdoh/blob/master/Report/Heatmap_Third.jpg)

From this step, 2 third order intervention is detected:

{'ACLS' } --> Upregulation, {'BTDDx'} --> Upregulation, {'PFOR'} --> Upregulation

{'ACLDC'} --> Upregulation,    {'ACLS' } --> Upregulation,  {'PFOR'} -->  Upregulation


One question is how other important metabolite fluxes change with these interventions. The following block answers that question:

First Order Intervention

PFOR Upregualtion:

```

uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [29];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [2];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)

>>>

T =

  7×18 table

                    Var1        Var2        Var3       Var4        Var5         Var6         Var7       Var8        Var9       Var10       Var11       Var12       Var13       Var14       Var15       Var16       Var17       Var18  
                  ________    ________    ________    _______    ________    __________    ________    _______    ________    ________    ________    ________    ________    ________    ________    ________    ________    ________

    EX_co          -20.042     -20.042     -20.042    -20.042     -20.042       -20.042     -20.042    -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042
    EX_co2          4.1653       3.846      2.4871     3.7134      3.3444       -1.0135      4.4383     4.2701       3.327      3.7676      3.2818      4.0642      3.6095      2.1347      4.4669      4.3336      3.6442      3.0801
    EX_h2          -33.042     -33.042     -33.042    -33.042     -33.042       -33.042     -33.042    -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042
    EX_ac           5.4925      4.3621      4.4716     4.2379      3.7611     0.0025362      3.8284     3.4852      4.3914      3.2449      5.1128       4.064      4.7102      5.3682      4.4285      4.4546      3.9932      5.3993
    EX_etoh         1.6866      2.4629      2.9538     2.2529      3.4321      0.069195      2.6112     2.9166      2.2006      3.7107      1.9565      2.5622      2.0058      2.1607      2.1265      2.2216      2.9281       2.028
    EX_bdoh        0.03979     0.10366    0.091711    0.15644    0.047319    1.1216e-07    0.081596    0.10593     0.10114     0.11823     0.14026    0.067976     0.14917     0.12378     0.16187     0.12737     0.10042     0.10046
    EX_biomass    0.017218    0.032191    0.029395    0.03154    0.022193    5.6673e-05    0.022987    0.02689    0.030683    0.016505    0.027154    0.026301    0.030764    0.033653    0.020412    0.025948    0.023239    0.027794


```

Second Order Interventions:

{'BTDDx'} --> Upregulation , {'PFOR' } --> Knockout

```
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [23,29];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [23,29];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)

>>>
T =

  7×18 table

                    Var1        Var2        Var3        Var4        Var5        Var6        Var7       Var8        Var9       Var10       Var11       Var12       Var13       Var14       Var15      Var16       Var17       Var18  
                  ________    ________    ________    ________    ________    ________    ________    _______    ________    ________    ________    ________    ________    ________    _______    ________    ________    ________

    EX_co          -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042    -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042    -20.042     -20.042     -20.042     -20.042
    EX_co2          4.1673      3.8596      2.6946      3.9309      3.3495      4.4229      4.5407     4.0297      3.3387      3.8498      3.4511      4.1745      3.6164      2.2294     4.4609      4.0128      3.6613      2.8454
    EX_h2          -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042    -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042    -33.042     -33.042     -33.042     -33.042
    EX_ac           5.5292      4.3538      4.3662      4.0911      3.7665      5.4039      3.8229     3.5861      4.3642      3.1911      4.8217      3.9774      4.7213      5.2922     4.4423      4.7759      3.9817      5.2464
    EX_etoh          1.649       2.457      2.9656      2.2932       3.423      1.3275      2.6175     2.9646      2.2233      3.7014      2.1256      2.5621        1.99      2.1797     2.1178      2.0761       2.934      2.2636
    EX_bdoh       0.038002     0.10394    0.088666      0.1644    0.047461    0.055243    0.081804    0.11179     0.10041     0.13082     0.13968    0.067929     0.14918     0.12362    0.16378     0.13302    0.099737     0.11924
    EX_biomass    0.017138    0.032465    0.028661    0.032146    0.022163    0.017299    0.020527    0.02632    0.030546    0.016832    0.026983    0.026439    0.031164    0.033669    0.02039    0.026331    0.023044    0.027845

```

{'ACLS'} --> Upregulation , {'PFOR' } --> Knockout

```

uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [59,29];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [2,2];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)


T =

  7×18 table

                    Var1        Var2        Var3        Var4        Var5         Var6         Var7        Var8        Var9       Var10      Var11       Var12       Var13       Var14       Var15       Var16       Var17      Var18  
                  ________    ________    ________    ________    ________    __________    ________    ________    ________    _______    ________    ________    ________    ________    ________    ________    _______    ________

    EX_co          -20.042     -20.042     -20.042     -20.042     -20.042       -20.042     -20.042     -20.042     -20.042    -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042    -20.042     -20.042
    EX_co2          4.1888      4.4295      2.4488      3.8516       3.419       -1.0407      4.7857       4.884      3.3689     3.7531      3.2824      4.1843      3.6113      2.1801      4.5689      4.5217     4.0806      3.1753
    EX_h2          -33.042     -33.042     -33.042     -33.042     -33.042       -33.042     -33.042     -33.042     -33.042    -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042    -33.042     -33.042
    EX_ac           5.4615       3.921      4.6959      4.1315      3.6915     0.0023367      3.6418      3.2444      4.2881     3.2299      5.1797      3.9948      4.6958      5.3761      4.3675      4.4277     3.7322      5.4313
    EX_etoh         1.6697       2.562      2.7915       2.244      3.4232       0.06547      2.6395       2.827      2.2235     3.6414      1.9044      2.5594      1.9328      2.0716       2.078      2.2175     2.9611      1.9177
    EX_bdoh       0.068341     0.16985     0.10084     0.17425    0.082228    2.0242e-07      0.1369     0.14045     0.12916    0.16478     0.16592    0.081936     0.19105     0.14283     0.16684     0.15559    0.12043     0.10139
    EX_biomass    0.016549    0.028017    0.027114    0.030284    0.020419    6.3542e-05    0.019464    0.023476    0.030482    0.01663    0.022738    0.026375    0.032335    0.032517    0.018991    0.025936    0.01833    0.025633


```

{'ACLS'} --> Upregulation , {'GLUD' } --> Knockout

```

uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [59,20];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [2,0];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)


T =

  7×18 table

                    Var1         Var2         Var3         Var4         Var5       Var6      Var7         Var8         Var9         Var10        Var11       Var12        Var13        Var14        Var15       Var16        Var17        Var18  
                  _________    _________    _________    _________    _________    ____    _________    _________    _________    _________    _________    ________    _________    _________    _________    ________    _________    _________

    EX_co           -20.042      -20.042      -20.042      -20.042      -20.042     0        -20.042      -20.042      -20.042      -20.042      -20.042     -20.042      -20.042      -20.042      -20.042     -20.042      -20.042      -20.042
    EX_co2           4.2927       5.7815       2.9838       3.1799       3.7625     0         4.6196       2.7956       4.0269       5.1427       3.4881      4.1912       4.6954       2.4943       6.0451      3.9626       5.4009       3.1758
    EX_h2           -33.042      -33.042      -33.042      -33.042      -33.042     0        -33.042      -33.042      -33.042      -33.042      -33.042     -33.042      -33.042      -33.042      -33.042     -33.042      -33.042      -33.042
    EX_ac            5.2558       3.3179       4.2396       4.6839       3.5368     0          4.195       4.1003       4.2821       2.9784        5.362      4.2683       3.8945       4.7528       3.9279      5.2861       3.4011       4.5474
    EX_etoh          1.9291       2.8292       3.3643        2.596       3.6886     0          2.476       3.3191       2.3648       3.4313       2.0532      2.5525       2.6876       2.9988       2.2103      1.8426       3.0226       3.0582
    EX_bdoh        0.059098      0.14503     0.048261      0.11055     0.065661     0        0.12375      0.17399     0.076534      0.13872     0.086496    0.083398      0.12445      0.13469       0.1297    0.044722     0.088716      0.13345
    EX_biomass    0.0076024    0.0063862    0.0071537    0.0088692    0.0087313     0      0.0039888    0.0079912    0.0062012    0.0017985    0.0055556    0.011488    0.0047478    0.0045674    0.0090837     0.01145    0.0023677    0.0062007
    
    
```

Third Order Interventions: 


{'ACLS'} --> Upregulation , {'BTDDx' } --> Upregulation, {'PFOR' } --> Upregulation

```
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [59,23,29];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [2,2,2];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)

T =

  7×18 table

                    Var1        Var2        Var3        Var4        Var5         Var6         Var7        Var8        Var9       Var10      Var11       Var12       Var13       Var14       Var15       Var16       Var17       Var18 
                  ________    ________    ________    ________    ________    __________    ________    ________    ________    _______    ________    ________    ________    ________    ________    ________    ________    _______

    EX_co          -20.042     -20.042     -20.042     -20.042     -20.042       -20.042     -20.042     -20.042     -20.042    -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042     -20.042    -20.042
    EX_co2          4.1888      4.4822      2.6966      4.1324      3.4245       -1.0177      4.8146       4.134      3.3886     3.8211      3.4769      4.2519      3.6338      2.2902      4.5667      4.7677      4.1335     2.6907
    EX_h2          -33.042     -33.042     -33.042     -33.042     -33.042       -33.042     -33.042     -33.042     -33.042    -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042     -33.042    -33.042
    EX_ac           5.5302      3.9091      4.5211      3.9398      3.6977     0.0023783      3.6696      3.5003      4.2533     3.1834      4.8729      3.9587      4.7055      5.2861      4.3785      4.2413      3.7098     5.1129
    EX_etoh          1.598      2.5402       2.844      2.2894      3.4127      0.068812      2.6284      2.8727      2.2518     3.6082      2.0772      2.5499      1.9092      2.0954      2.0706       2.269      2.9665     2.4028
    EX_bdoh       0.066847      0.1752    0.097231     0.18372    0.082336     2.169e-07     0.13787     0.14263     0.12878    0.19042     0.16641    0.081614     0.19794     0.14349     0.16917     0.15313     0.11952    0.16626
    EX_biomass    0.016744    0.028076    0.026014    0.030976    0.020362    8.5448e-05    0.018303    0.023115    0.030262    0.01676    0.022452    0.026482    0.033195    0.032687    0.018983    0.026182    0.017607    0.02557


```


{'ACLDC'} --> Upregulation , {'ACLS' } --> Upregulation, {'PFOR' } --> Upregulation

```
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};Uptakes.values{2} ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

perturbed_rxn = [13,59,29];                                                         % rxn # for enzyme to be changed; = [] if no enzyme change
expression_level = [2,2,2];                                                      % = 1 if no enzyme change; = 0 if knockout; > 1 if overexpression; < 1 if underexpression
t_interval = 0:5e2;                                                         % INPUT - time interval for DAE integration

    save_results = 0;                                                           % 1-save results as matlab file; 0-do not save results
    if save_results
        file_save_name = 'TEST_results';                                        % INPUT desired name of results file
    end



    
    [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
        perturb_Ksets(model, All_K_final, All_fractions_final, ...
        t_interval, uptake_rxns, uptake_values, perturbed_rxn, ...
        expression_level);

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)

T =

  7×18 table

                    Var1        Var2        Var3        Var4        Var5         Var6         Var7        Var8        Var9       Var10       Var11       Var12      Var13       Var14       Var15       Var16      Var17       Var18  
                  ________    ________    ________    ________    ________    __________    ________    ________    ________    ________    ________    _______    ________    ________    ________    _______    ________    ________

    EX_co          -20.042     -20.042     -20.042     -20.042     -20.042       -20.042     -20.042     -20.042     -20.042     -20.042     -20.042    -20.042     -20.042     -20.042     -20.042    -20.042     -20.042     -20.042
    EX_co2          4.1905      4.4478      2.4481      3.7591      3.4272       -1.0106      4.8214      2.3621      3.3676      3.7637      3.2953     4.3953      3.6672      2.0866      4.6694     4.6177      4.1981        2.84
    EX_h2          -33.042     -33.042     -33.042     -33.042     -33.042       -33.042     -33.042     -33.042     -33.042     -33.042     -33.042    -33.042     -33.042     -33.042     -33.042    -33.042     -33.042     -33.042
    EX_ac           5.4799      3.8991      4.6969      4.2391      3.6994     0.0026326      3.6301      3.5657      4.3056       3.212      5.2087     3.7897      4.6751      5.3737      4.3788     4.3167      3.7062      5.3384
    EX_etoh         1.6534      2.5263       2.791      2.1893      3.4209      0.070363      2.6365      2.5398      2.2033      3.5871      1.8856     2.5813      1.9139      2.0398      2.0315     2.2442      2.9664      2.1332
    EX_bdoh       0.065137     0.21521      0.1017     0.22474    0.082684    3.5102e-07     0.14371     0.14774     0.13209     0.20585     0.17924    0.09303     0.22044     0.22712     0.24381    0.16407     0.14966     0.14736
    EX_biomass    0.016099    0.028076    0.026976    0.026211    0.020111    5.1872e-05    0.018873    0.024167    0.029684    0.017615    0.022007    0.02574    0.030036    0.031587    0.019522    0.02621    0.018714    0.025064


