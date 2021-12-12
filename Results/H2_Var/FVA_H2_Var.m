Metaclau=readCbModel('Metaclau.mat')
load('EXP_Struct.mat')

for i=1:size(EXP_Struct,2)
    Uptakes.rxns{i}=EXP_Struct(i).Name;
    Uptakes.values{i}=EXP_Struct(i).EXP1;                                     %If you want another experiment, add that to EXP_Struct, and change EXP1 accordingly in this line
    Uptakes.type{i}='b';
end
Metaclau=changeRxnBounds(Metaclau,Metaclau.rxns(789),Uptakes.values{1});
Metaclau.rxns(find(Metaclau.c))
inds=[788,771,772,790,786,787,849]
i=1;
for Hyd_Ups=0:-0.5:-50
    Temp_Meta=changeRxnBounds(Metaclau,'EX_HYDROGEN-MOLECULE',Hyd_Ups,'l')
    [min,max]=fluxVariability(Temp_Meta, 'optPercentage', 75, 'allowLoops',1,'printLevel',2,'rxnNameList',Temp_Meta.rxns(inds),'method','1-norm');
    Min_Sol(:,i)=min;
    Max_Sol(:,i)=max;
    i=i+1;

end

for J=1:size(Min_Sol,1)
    subplot(size(Min_Sol,1),1,J)
    plot(-0:-0.5:-50,Min_Sol(J,:),'.',0:-0.5:-50,Max_Sol(J,:),'.')
    legend({'Min','Max'})
    title(regexprep(Metaclau.rxns(inds(J)),'_','-'))
    xlabel('Hydrogen Uptake Flux')
    ylabel('Exchange Flux')
end

