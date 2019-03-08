function myfitness=runODEmodel(bigvaccine,problem_info, use_all_antigens, objective, use_weighted, use_one_antigen, returnfull)
% function myfitness=getInvAny(bigvaccine,problem_info, use_all_antigens, objective, use_weighted, use_one_antigen, returnfull)
% this function computes the invasivness or DR-inv combo under flexible options: 
% (1) vary all antigens or not (1 or 0 in use_all_antigens)
% (2) objective: which objective is computed. If 'returnfull' is 0, only this number is returned. Specify 'invexp',
% 'drexp' or 'kidsexp' for invasiveness (combo), combined inv and dr or kid only under the exponential model. 
% (2) weighted or not weighted (different strenghts of NFDS)
% (3) a SINGLE antigen specified by an integer in 1 to 12 (MA) or 1:10 (Maela) (use_one_antigen=4 to use the 4th one, etc)
% (4) If 'returnfull' is 1, return the solution to the ODE and some other information including all the objectives. 
% 'returnfull' is set to 0 to use this in optimisation routines.
% ODE options are hard coded in this function
% 10 year final time frame is hard coded too
% the 'vacfactor' of 0.06 is hard coded (it was v=0.1 in Nick's model; this fit test conditions best to match Nick's model)

%% set up parameters
numSeros=length(problem_info.seronames); % 32 or 63 (maela), last one is 'NT' 
LL= size(problem_info.G,1); % number of isolates
sigma=problem_info.sigma;
ops=odeset('AbsTol',1e-8,'RelTol',1e-5,'NonNegative', 1:LL,'BDF','off');
tmax=12*10; %% NOTE HARD CODED that we solve the model for ONLY 10 YEARS. 

% antifactor=0.05; % how well do the antigens do, compared to the serotypes? st have 0.1. now half that. 

if use_weighted==1
    % vacfactor=0.0812; 
    vacfactor=problem_info.v; % was 0.0812 as above
    antifactor=0.5*problem_info.v;
else
    vacfactor=0.06; % was 0.1 in Nick's model but my fits suggested 0.06; anyway I use the weighted case now
end
m=problem_info.m;


%% set up serotype part of the vaccine strategy 
serovt=bigvaccine(1:(numSeros-1));
serovt(numSeros)=0; % don't vaccinate the NT (non-typables)
VT=problem_info.SerotypeToStrain*serovt';  % serovt: is serotype j in the vaccine?
VT(VT>0)=1; 

if use_all_antigens==0
    combiVT=vacfactor*VT;
end


%% set up antigen part of the vaccine strategy 
if use_all_antigens==1
   antigenvaccine=bigvaccine(numSeros:end); 
   antigenvals=problem_info.antigenvals; 
   antigenvals(problem_info.AntiINDEX(antigenvaccine==0))=0; % no antigen effect where the antigen is not in the vaacine    
   antigeneffect = max( problem_info.G .* repmat(antigenvals, LL, 1), [], 2); % max of all antigen effects in each isolate
   combiVT=max(vacfactor*VT, antifactor*antigeneffect);
end

if use_all_antigens==1 & use_one_antigen==1
    error('use either all antigens or only one, but not both')
end

if use_one_antigen > 0 % then its value says which one of the 12 (or 10 maela) is included
   
   antigenvals=problem_info.antigenvals; 
   antigenvaccine=zeros(1,length(problem_info.AntiINDEX));
   antigenvaccine(use_one_antigen)=1; 
   % set all but the ONE in use_one_antigen to 0  
   antigenvals(problem_info.AntiINDEX(antigenvaccine==0))=0; % no antigen effect where the antigen is not in the vaacine    
   thislocus=problem_info.AntiINDEX(use_one_antigen); 
   antigeneffect=problem_info.G(:,thislocus)*problem_info.antigenvals(thislocus); 
   combiVT=max(vacfactor*VT, antifactor*antigeneffect);
end


%% set up ode function handle and initial conditions
rr=-log(1-combiVT); % from first order discretization (in reverse): combiVT includes Nick/Jukka's 'v' parameter 
rho=log(1+sigma);  % from first order discretization (in reverse); sigma is comparable to Nick/Jukka choice of sigma
K=100000; % carrying capacity

% use specified initial conditions if available
if any(strcmp(fieldnames(problem_info),'ics'))
    x0=problem_info.ics;
else
    x0(VT==1)=0.5; x0(VT==0)=0.5; x0=x0*(K/sum(x0));
    % everything is equally prevalent at the start if we don't pass in any information
end

if use_weighted==0
     diffxt=@(t,x) x ...
        .*(log(K/sum(x)) - rr + rho*(problem_info.dev-problem_info.G*problem_info.G'*(x/sum(x)))) ;
    
else % version for variable strength of nfds
    modifiedGT=repmat(problem_info.locusweights, 1, size(problem_info.G,1)).*problem_info.G';
% this ODE comes from the standard discretization (in reverse)
    diffxt=@(t,x) x ...
        .*(log(K/sum(x)) - rr + rho* (problem_info.weighteddev-problem_info.G*modifiedGT*(x/sum(x)))) + m*ones(size(x)) ;
end
  
%% solve ode and compute fitness cost 
[t,xsol]=ode15s(diffxt,[0 tmax], x0,ops);  % def not at eq after 10000 mathematically but invasiveness equilibrates soon
myfitness = getInvasiveness(xsol(end,:), problem_info, objective);

%% create detailed output if requred
if returnfull==1
oderesult.t=t;
oderesult.xsol=xsol;
oderesult.fitness=myfitness;
oderesult.VT=VT; 

oderesult.Inv=getInvasiveness(xsol(end,:), problem_info,'Inv');

oderesult.InvKids=getInvasiveness(xsol(end,:), problem_info,'KidsOnly');

oderesult.InvAdults=getInvasiveness(xsol(end,:), problem_info,'AdultsOnly');

oderesult.InvWithDR=getInvasiveness(xsol(end,:), problem_info,'DR');

oderesult.kidsexp=getInvasiveness(xsol(end,:), problem_info,'kidsexp');
oderesult.invexp=getInvasiveness(xsol(end,:), problem_info,'invexp');
oderesult.drexp=getInvasiveness(xsol(end,:), problem_info,'drexp');

oderesult.DR=xsol(end,:)*problem_info.DRscore/sum(xsol(end,:));
oderesult.scfreq=groupbySC(xsol(end,:), problem_info); 
oderesult.serovt=serovt; 
if use_all_antigens==1; oderesult.antigenvaccine=antigenvaccine; end 
if use_one_antigen >0 ; oderesult.use_one_antigen=use_one_antigen; end
myfitness=oderesult;
end

