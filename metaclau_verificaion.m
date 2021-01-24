load('Metaclau.mat')
metaclau.rxns
Rxns=[{'EX_CARBON-MONOXIDE'},
      {'EX_HYDROGEN-MOLECULE'},
      {'EX_ETOH'},                                       
      {'EX_BUTANEDIOL'}];
 
inds=find(ismember(metaclau.rxns,Rxns));
metaclau=changeRxnBounds(metaclau,'EX_CARBON-MONOXIDE',-16.5,'l')
Res=optimizeCbModel(metaclau,'max','one')
metaclau=changeRxnBounds(metaclau,'EX_BIOMASS',Res.f);
metaclau=changeRxnBounds(metaclau,'EX_CARBON-MONOXIDE',-18.4,'b')
optimizeCbModel(metaclau,'max','one')