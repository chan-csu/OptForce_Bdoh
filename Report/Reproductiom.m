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


%% bdoh maximizer mutant
MT=changeObjective(metaclau,'EX_BUTANEDIOL')
solution=optimizeCbModel(MT,'max','one');
MT_exports=solution.x(metaclau_Indices)
table(metaclau.rxns(metaclau_Indices),MT_exports)



