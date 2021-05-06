Metaclau=readCbModel('Metaclau.mat')
iCLAU=readCbModel('iCLAU786_2018.xls')
% Medium=

%% Test 1 Flux Through the network with zero Exchanges
Meta_Exc=findExcRxns(Metaclau);
Meta_Uptakes=find(Metaclau.lb(Meta_Exc)<0);
Metaclau_Test1=changeRxnBounds(Metaclau,Metaclau.rxns(Meta_Exc),0,'b');

i_Exc=findExcRxns(iCLAU);
i_Uptakes=find(iCLAU.lb(i_Exc)<0);
iCLAU_Test1=changeRxnBounds(iCLAU,iCLAU.rxns(i_Exc),0,'b');
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






%% Test 3 Matching the experimental data