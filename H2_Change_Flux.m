load Final_Kinetic_Ensemble.mat
load('model_core.mat')
for i=[-5 -10 -15 -20 -25 -30 -35 -40 -45 -50 -55 -60 -65 -70 -75 -80 -85 -90 -95 -100]
uptake_rxns = [1;2]; %[1 - CO Uptake; 2 - H2 Uptake]                        % INPUT elementary rxn index for uptake rxns constrained for DAE 
uptake_values = [Uptakes.values{1};i ];                                           % INPUT uptake fluxes for CO and H2 gases.  Current values are experimental values at high biomass concentration.

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

T=array2table(solutions);
T.Properties.RowNames=model.rxns(:)
writetable(T,'Kinetic_Model_Fluxes_CO_Constant.xlsx','Sheet',num2str(i),'WriteRowNames',true)
end