%% Metaclau Wild Type
load('EXP_Struct.mat')

for i=1:size(EXP_Struct,2)
    Uptakes.rxns{i}=EXP_Struct(i).Name;
    Uptakes.values{i}=EXP_Struct(i).EXP1;                                     %If you want another experiment, add that to EXP_Struct, and change EXP1 accordingly in this line
    Uptakes.type{i}='b';
    metaclau=changeRxnBounds(metaclau,Uptakes.rxns{i},Uptakes.values{i},Uptakes.type{i})
end



metaclau_Indices=[789,788,790,771,786,787,849];

Biomass_ID=849;

Solution=optimizeCbModel(metaclau,'max','one');
WT_exports=Solution.x(metaclau_Indices)
table(metaclau.rxns(metaclau_Indices),WT_exports)

%% Kinetic Model Wild type

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

Kin_Model_Indices=[1 3 2 58 54 62 64 ]
T=array2table(solutions(Kin_Model_Indices,:))
T.Properties.RowNames=model.rxns(Kin_Model_Indices)
Base_BDOH=solutions(62,:);

%% bdoh maximizer mutant
MT=changeObjective(metaclau,'EX_BUTANEDIOL')
solution=optimizeCbModel(MT,'max','one');
MT_exports=solution.x(metaclau_Indices)
table(metaclau.rxns(metaclau_Indices),MT_exports)

%% Stat Analysis of the interventions

%selection of first order OptForce rxns that are in core model
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

[p,tbl,stats] = anova1(First_order_results(:,:,3));
[c,m,h,gnames] = multcompare(stats);


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
