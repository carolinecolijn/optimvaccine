function result=getInvErrorbars(xsolend, dataset, nReps)
% function InvSample=getRandomInvasiveness(xsolend, dataset,objective)
% idea: for each isolate: sample inv for kids , sample inv for adults, take mean. 

% samples are normal with appropriate CIs as in Nick's meta-analysis
% doesn't allow 'kids only', 'adults only' if we take the mean; may want to switch on an 
% objective, kids, adults, etc. and maybe have a switch to use randomness or not? 
if length(xsolend)~=size(dataset.G,1)
    error('trying to compute error based on a simulation that does not match the input dataset')
end

% plan: 

% reprogram efficiently if this gets too slow, but I doubt it'll be the limiting factor

% pull st-specific values of the objectives and their sigmas from the genotype-expanded ones in the dataset: 
for n=1:length(dataset.seronames)
    myinds=find(strcmp(dataset.serotype, dataset.seronames{n})); % indices of genomes with this ST
    adultinv(n)=dataset.Invasiveness.adults(myinds(1));
    kidinv(n)=dataset.Invasiveness.kids(myinds(1));
    adultsig(n)=dataset.Invasiveness.adultsSigma(myinds(1)); % value for that serotype in solution X
    kidssig(n)=dataset.Invasiveness.kidsSigma(myinds(1)); % value for that serotype in solution X
end


 nvec=xsolend(:);
 % set up for DRexp objective
 a= -2; b= 0.5;  myprobs = (1./ (1+exp(-a-b*dataset.DRscore)));
totalPop=sum(nvec); 
 newadultgenoinv=zeros(size(nvec));
 newkidgenoinv=zeros(size(nvec)); 
 % do the following many times:
 for k=1:nReps;
     
     % --- sample a new invasiveness for each serotype
     newadultseroinv=normrnd(adultinv, adultsig);
     newkidseroinv=normrnd(kidinv, kidssig);
     % --- give each individual in the population that invasiveness
     for n=1:length(dataset.seronames)
         myinds=find(strcmp(dataset.serotype, dataset.seronames{n}));
         newadultgenoinv(myinds)=newadultseroinv(n);
         newkidgenoinv(myinds)=newkidseroinv(n);
     end
     
     % --- add up the total invasiveness in the population
     totcombo(k)=sum(nvec(:).*exp(0.5*(newadultgenoinv(:)+newkidgenoinv(:))))/totalPop;
     totkids(k)=sum(nvec.*exp(newkidgenoinv))/totalPop;
     newObj = exp(0.5*(newadultgenoinv+newkidgenoinv)).*myprobs; % each genotype's contribution to the total
     totdr(k) = sum( nvec(:).*newObj(:))/totalPop;
     
 end
 % then look at lower and upper quartiles for ex.

        
 result.ranges=[ quantile(totdr,0.25) quantile(totdr,0.75) ...
     quantile(totcombo,0.25) quantile(totcombo,0.75) ...
     quantile(totkids,0.25) quantile(totkids,0.75)];
     
result.totcombo=totcombo;
result.totkids=totkids;
result.totdr=totdr; 
result.note='dr 0.25, dr 0.75, inv 0.25, inv 0.75, kids 0.25, kids 0.75'; 

% IDEA: we could re-run EVERYTHING using as an objective the probability that 
% the invasivness does better than (whatever). 
% or otherwise incorporate this better into the optimisation. 
% oh, the things you learn. 





