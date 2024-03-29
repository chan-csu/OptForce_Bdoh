% 
% initCobraToolbox(false)
% 
% solverStatus = changeCobraSolver('gurobi', 'all'); 

Metaclau=readCbModel('Metaclau.mat');
iCLAU=readCbModel('iCLAU786_2018.xls');

Metaclau=changeObjective(Metaclau,Metaclau.rxns(849));
iCLAU=changeObjective(iCLAU,iCLAU.rxns(989));

% Medium=

%% Test 1 Flux Through the network with zero Exchanges
clc
fprintf('-------***Test 1***-----------------\n')
Meta_Exc=findExcRxns(Metaclau);
Meta_Uptakes=find(Metaclau.lb(Meta_Exc)<0);
Metaclau_Test1=changeRxnBounds(Metaclau,Metaclau.rxns(Meta_Exc),0,'b');
Metaclau_Test1.c=ones(size(Metaclau_Test1.rxns));
i_Exc=findExcRxns(iCLAU);
i_Uptakes=find(iCLAU.lb(i_Exc)<0);
iCLAU_Test1=changeRxnBounds(iCLAU,iCLAU.rxns(i_Exc),0,'b');
iCLAU_Test1.c=ones(size(iCLAU_Test1.rxns));
Sol_Meta_T1=optimizeCbModel(Metaclau_Test1);
Sol_i_T1=optimizeCbModel(iCLAU_Test1);

if isempty(Sol_Meta_T1.x~=0)
    disp('Metaclau passed Test 1')
    
else
    disp('Metaclau failed Test 1: The following reactions had non-zero fluxes:')
    fprintf('%s',Metaclau.rxns(find(Sol_Meta_T1.x~=0)))
end


if isempty(Sol_i_T1.x~=0)
    disp('iCLAU passed Test 1')
    
else
    disp('iCLAU failed Test 1: The following reactions had non-zero fluxes:')
    fprintf('%s',iCLAU.rxnNames{find(Sol_i_T1.x~=0)})
end

%% Test 2 Flux Through ATPm with zero uptakes
fprintf('\n-------***Test 2***-----------------\n')
Metaclau_Test2=changeObjective(Metaclau_Test1,'ATPASE-RXN');
Sol_Meta_T2=optimizeCbModel(Metaclau_Test2);

if Sol_Meta_T2.f==0 | Sol_Meta_T2.origStat=='INFEASIBLE'
    fprintf('\nMetaclau passed Test 2')
    
else
    fprintf('\nMetaclau failed Test 2: Non-zero flux through ATPm')
end

iCLAU_Test2=changeObjective(iCLAU_Test1,'rxn00062_c0');
Sol_i_T2=optimizeCbModel(iCLAU_Test2);
if Sol_i_T2.f==0 | strcmp('INFEASIBLE',Sol_i_T2.origStat)
    fprintf('\niclau passed Test 2')
    
else
    fprintf('\niCLAU failed Test 2: Non-zero flux through ATPm\n')
end


%% Test 3 Predicting the experimental data
%Here we test against three different Experimental datasets
EXP=readtable('Exp_Data.csv','TreatAsEmpty',{'NA'});
%the following indicates indexing for exchange reaction in both models. The
%order is the same as column order in 

Metaclau_Inds=[779,782,781,784,780,798,789,788,790,771,786,787,849];
iCLAU_Inds=[61,1008,852,100000,51,406,960,854,966,972,969,986,989];
Metaclau_Test3=changeRxnBounds(Metaclau,Metaclau.rxns(Meta_Uptakes),-1,'l');
Metaclau_Test3=changeRxnBounds(Metaclau_Test3,Metaclau.rxns(Metaclau_Inds(1:6)),-1000,'l');
iCLAU_Test3=changeRxnBounds(iCLAU,iCLAU.rxns(i_Uptakes),-1,'l');
iCLAU_Test3=changeRxnBounds(iCLAU_Test3,iCLAU.rxns(iCLAU_Inds([1 2 3 5 6])),-1000,'l');
NaN_Inds=find(any(isnan(table2array(EXP)), 2))';
for i=1:17
    Temp_Meta=changeRxnBounds(Metaclau_Test3,Metaclau_Test3.rxns(Metaclau_Inds([7:9])),EXP{i,[1:3]},'b');
    tmp=optimizeCbModel(Temp_Meta,'max','one');  
    if tmp.f == 0 | strcmp('INFEASIBLE',tmp.origStat)
        Sim_Res_Meta(i,:)= [0 0 0 0];
        NaN_Inds=[NaN_Inds i];
    else
        Sol_Meta_3(i)=optimizeCbModel(Temp_Meta,'max','one');
        Sim_Res_Meta(i,:)=Sol_Meta_3(i).x(Metaclau_Inds(1,10:13));
    end
    
    Temp_i=changeRxnBounds(iCLAU_Test3,iCLAU_Test3.rxns(iCLAU_Inds([7:9])),EXP{i,[1:3]},'b');
    tmp=optimizeCbModel(Temp_i,'max','one');
    if tmp.f==0 | strcmp('INFEASIBLE',tmp.origStat)
        Sim_Res_i(i,:)= [NaN NaN NaN NaN];
        NaN_Inds=[NaN_Inds i];
    else
        Sol_i_3(i)=optimizeCbModel(Temp_i,'max','one');
        Sim_Res_i(i,:)=Sol_i_3(i).x(iCLAU_Inds(10:13)');
    end
end
Sim_Res_Meta(NaN_Inds,:)=[];
Sim_Res_i(NaN_Inds,:)=[];
EXP(NaN_Inds,:)=[];
fprintf('\n-------***Test 3***-----------------\n')
fprintf('The norm of difference between Predicted data and experimental data for Metaclau is: %s \n ',norm(Sim_Res_Meta-table2array(EXP(:,4:end))))
fprintf('\nThe norm of difference between Predicted data and experimental data for iCLAU is: %s ',norm(Sim_Res_i-table2array(EXP(:,4:end))))


for i=1:size(EXP,1)
    fprintf('\n----***FVA results for Metaclau***----\n')
    Temp_Meta=changeRxnBounds(Metaclau_Test3,Metaclau_Test3.rxns(Metaclau_Inds([7:9])),EXP{i,[1:3]},'b');
    Sol = optimizeCbModel(Temp_Meta,'max','one');
    if Sol.f==0 | strcmp('INFEASIBLE',Sol.origStat)
        Sim_Res_Meta_min(i,:) = [NaN NaN NaN NaN];
        Sim_Res_Meta_max(i,:) = [NaN NaN NaN NaN];
    else
        [fluxMin, fluxMax] = fluxVariability(Temp_Meta, 'optPercentage', 75, 'allowLoops',1,'printLevel',2,'rxnNameList',Metaclau_Test3.rxns(Metaclau_Inds(1,10:13)),'method','1-norm');
        rxns_min_model_check = fluxMin;
        rxns_max_model_check = fluxMax;
        Sim_Res_Meta_min(i,:) = rxns_min_model_check;
        Sim_Res_Meta_max(i,:) = rxns_max_model_check;
    end
    
    fprintf('\n----***FVA results for iCLAU***----\n')
    Temp_i=changeRxnBounds(iCLAU_Test3,iCLAU_Test3.rxns(iCLAU_Inds([7:9])),EXP{i,[1:3]},'b');
    Sol = optimizeCbModel(Temp_i,'max','one');
    if Sol.f==0 | strcmp('INFEASIBLE',Sol.origStat)
        Sim_Res_Meta_min(i,:) = [NaN NaN NaN NaN];
        Sim_Res_Meta_max(i,:) = [NaN NaN NaN NaN];
    else
    [fluxMin, fluxMax] = fluxVariability(Temp_i, 'optPercentage', 75, 'allowLoops',1,'printLevel',2,'rxnNameList',Temp_i.rxns(iCLAU_Inds(1,10:13)),'method','1-norm');
    rxns_min_model_check = fluxMin;
    rxns_max_model_check = fluxMax;
    Sim_Res_i_min(i,:) = rxns_min_model_check;
    Sim_Res_i_max(i,:) = rxns_max_model_check;
    end
end



% R_cpd00009_ext_b Phosphate exchange
% R_cpd00013_ext_b   Ammonia Exchange
% R_cpd00001_ext_b Water                
% R_cpd00011_ext_b CO2 Exchange
% cpd00067_ext_b Proton Exchange
% R_cpd00204_ext_b CO Exchange
% R_cpd11640_ext_b H2 Exchange
% R_cpd00363_ext_b Ethanol Exchange
% R_cpd00029_ext_b Acetate exchange
% R_cpd01947_ext_b BDOH
% R_cpd11416_ext_b Biomass
% cpd00239_ext_b for H2S
%     
% No NIACIN 100000 as index
