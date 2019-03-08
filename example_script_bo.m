

Basic={'1',  '5', '14'} ; % ALWAYS
nAllSeros=length(massdata.seronames)-1; 

% set parameters for the model
use_weighted=1; 
use_one_antigen=0;
use_all_antigens=0; 
thisObjective = 'kidsexp'; 

% set up the serotypes
includedSTs=Basic; % high-invasiveness STs 
excludedSTs=[];
varyingSTs=setdiff(massdata.seronames(1:nAllSeros), ...
    union(includedSTs,excludedSTs))'; 

% set parameters, constraints and function handles for the BO
MaxHours=23;
dataset='mass';
MaxSerotypes=12; % including 1, 5 and 14 this means we would have 15-valent formulations


myobjectivefun= @(x) fitness_for_bo(x, includedSTs,excludedSTs, varyingSTs, massdata, use_all_antigens,  ...
    thisObjective, use_weighted, use_one_antigen); % objective function; solves the ODE

myxfun= @(x) myxconstraint(x, MaxSerotypes);  % constraint function

% create the variables matlab needs ('optimizableVariable' function) with integer type
NumAntigens=length(massdata.AntiINDEX); 

numVars=length(varyingSTs)+ NumAntigens*use_all_antigens; 
allvars=makeallvars(numVars); 

boresult=bayesopt(myobjectivefun,allvars,'XConstraintFcn',myxfun, ...
   'MaxObjectiveEvaluations',8000, 'MaxTime',0.001*MaxHours*60*60,'Verbose',2)

