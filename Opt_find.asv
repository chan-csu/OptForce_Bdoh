function [optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets,Eff_out]=Opt_find(Model,Uptakes)
if nargin==0
    Model=readCbModel('Metaclau.mat');
end
% Uptakes is a struct. Uptakes.rxns is a cell containing the names of the reactions we want to constrain
%Uptakes.values is the constrained flux of the reactions. it's a matrix
%containing numbers. Uptakes.types is the type of constraints. it is a cell
%contanining characters 'l', 'b', and 'u'
%% Setting up the Wild_Type and Mutant Solns
Model=readCbModel('Metaclau.mat');
for i=1:length(Uptakes.rxns)
    Model=changeRxnBounds(Model,Uptakes.rxns{i},Uptakes.values{i},Uptakes.type{i});
end
Model=changeObjective(Model,'EX_BIOMASS');
%   This can be changed to a conditional form: If solver is not already
%   gurobi, change it to gurobi
changeCobraSolver('ibm_cplex','all');
Wild_Type=optimizeCbModel(Model,'max');
Model_MT=changeObjective(Model,'EX_BUTANEDIOL');
Model_MT=changeRxnBounds(Model_MT,'EX_BIOMASS',Wild_Type.f*0.1,'b');
Mutant=optimizeCbModel(Model_MT);
exchanges=cellfun(@isempty, strfind(Model.rxns, 'EX_'))==0;
Ex_Rxns=Model.rxns(exchanges);Ex_MT=Mutant.x(exchanges);Ex_WT=Wild_Type.x(exchanges);
Diff_Ind=Ex_MT~=Ex_WT;
Exchange_Comp=table(Ex_Rxns(Diff_Ind),Ex_WT(Diff_Ind),Ex_MT(Diff_Ind),'VariableNames',{'Ex_Rxns','WT','MT'})
clear Ex_MT Ex_WT Diff_Ind

%% Building The Constraints
constrWT = struct('rxnList', {{'EX_BIOMASS'}}, 'rxnValues', Wild_Type.f, 'rxnBoundType', 'b');
constrMT = struct('rxnList', {{'EX_BIOMASS', 'EX_BUTANEDIOL'}}, 'rxnValues', [Wild_Type.f*0.1-0.00001, Mutant.f-0.000001],'rxnBoundType', 'bb');
[minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(Model,constrWT, constrMT);
FVA_Results=table(Model.rxns,minFluxesW, maxFluxesW, minFluxesM, maxFluxesM)
%% Finding Musts
runID = 'OptForce_bdoh';
constrOpt = struct('rxnList', {{'EX_BIOMASS','EX_BUTANEDIOL'}}, 'values', [Wild_Type.f*0.1-0.00001, Mutant.f-0.000001]');
[mustLSet, pos_mustL] = findMustL(Model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
                                  'outputFileName', 'MustL' , 'printExcel', 1, 'printText', 0, ...
                                  'printReport', 0, 'keepInputs', 1, 'verbose', 0);
[mustUSet, pos_mustU] = findMustU(Model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
                                  'outputFileName', 'MustU' , 'printExcel', 0, 'printText', 0, ...
                                  'printReport', 1, 'keepInputs', 1, 'verbose', 0);
                              
exchangeRxns = Model.rxns(cellfun(@isempty, strfind(Model.rxns, 'EX_')) == 0);
excludedRxns = unique([mustUSet; mustLSet; exchangeRxns])

[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
    findMustUU(Model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
               'printExcel', 1, 'printText', 0, 'printReport', 0, 'keepInputs', 1, ...
               'verbose', 1);
[mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
    findMustLL(Model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
               'printExcel', 1, 'printText', 0, 'printReport', 0, 'keepInputs', 1, ...
               'verbose', 1);
[mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = ...
    findMustUL(Model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUL', 'outputFileName', 'MustUL', ...
               'printExcel', 1, 'printText', 0, 'printReport', 0, 'keepInputs', 1, ...
               'verbose', 1);

mustUSet
mustLSet
mustUU
mustUL
mustLL
load 'Mapping vect.mat'
%% OptForce 
mustU = unique(union(mustUSet, mustUU));
mustL = unique(union(mustLSet, mustLL));
mustU=mustU(ismember(mustU,Mapping_Vect(:,1)));
mustL=mustL(ismember(mustL,Mapping_Vect(:,1)));
targetRxn = 'EX_BUTANEDIOL';
biomassRxn = 'EX_BIOMASS';
constrOpt = struct('rxnList', {{'EX_BIOMASS'}}, 'values', [0.1*Wild_Type.f]);
for i=1:3
[optForceSets{i}, posOptForceSets{i}, typeRegOptForceSets{i}, flux_optForceSets{i}] = ...
    optForce(Model, targetRxn, biomassRxn, mustU, mustL, minFluxesW, maxFluxesW, minFluxesM, maxFluxesM,'k',i,'nSets',500,'constrOpt', constrOpt);
end
%% Showing the final results
for i=1:3
    fprintf('Number of interventions = %d \n',i)
table(optForceSets{i}, posOptForceSets{i}, typeRegOptForceSets{i}, flux_optForceSets{i},'VariableNames',{'optForceSets' ...
    'indices','Type_of_Intervention','Fluxes'})
end
 Duble_int=posOptForceSets{2};
 j=0;
for i=1:size(Duble_int)
    if length(Model.rules{Duble_int(i,1)})~=0 && length(Model.rules{Duble_int(i,2)})~=0
        j=j+1;
        Eff_out(j,:)=Duble_int(i,:);
    end
end
end

