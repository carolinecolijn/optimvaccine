function objective = fitness_for_bo(x, includedSTs, excludedSTs, varyingSTs, problem_info, use_all_antigens, objective, use_weighted, use_one_antigen)
% function objective = fitness_for_bo(x, includedSTs, excludedSTs, varyingSTs, problem_info, use_all_antigens, objective, use_weighted, use_one_antigen)
% This function takes in x (in the format required by bayesopt) and parameters connecting x to the serotype names and data
% structure, and then computes the specified fitness of the bacterial population at the 10 year time point. 


% convert x to an array. 
arrx=table2array(x); 

% check that dimensions work. x has one component for each varying ST or antigen
if size(arrx,2) ~= length(varyingSTs) + use_all_antigens*length(problem_info.AntiINDEX)
    error('mismatch size between the optimisation variables and the fixed serotypes and antigens')
end 

% need to create a vector 'bigvaccine'. it has length (number of STs in the model that are not 'NT') + number of antigenic proteins
% so for the MA data it must be length 31, or 31+12 = 43 if we are using all 12 antigens. 
nAllSeros=length(problem_info.seronames)-1; % last is 'NT' in both MA and Maela data
serovt=zeros(1,nAllSeros);
 
% set 1 for all included STs
for n=1:length(includedSTs)
     serovt(strcmp(problem_info.seronames,includedSTs{n}))=1;
end

% set 0 for all excluded STs
for n=1:length(excludedSTs)
     serovt(strcmp(problem_info.seronames,excludedSTs{n}))=0;
end

% set according to x for the varying STs
for n=1:length(varyingSTs)
     serovt(strcmp(problem_info.seronames,varyingSTs{n}))=arrx(n); 
end

bigvaccine=serovt; 

% set antigen part according to x
if use_all_antigens==1
    antipart=arrx((length(varyingSTs)+1) : end );
    bigvaccine=[bigvaccine antipart]; 
end 


returnfull=0;

% run the ODE
objective=runODEmodel(bigvaccine,problem_info, use_all_antigens, objective, use_weighted, use_one_antigen, returnfull);













