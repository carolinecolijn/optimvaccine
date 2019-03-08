function InvSample=getInvasiveness(xsolend, problem_all, objective)
% function InvSample=getRandomInvasiveness(xsolend, problem_all,objective)
% inputs are the end point (or some time point) of the solution to the ODE, namely estimated
% prevalence of each unique genotype at a specific time, the data structure for hte problem, and the name of the
% objective. These can be: 
% (1) for kis invasiveness, use 'Exp','exp','kidsexp','KidsExp'
% (2) for combined adults and kids use 'InvExp','invexp','Invexp','invExp'
% (3) for dr, adults and kids use 'DRExp', 'drexp','Drexp','DRexp'
% this code doesn't have my approach for random invasivness. That's in getInvErrorbars.m
% parameters in the way we combined the DR score and the invasivness are hard coded here.         
        
 nvec=xsolend(:);
 
switch objective
   
      case {'Exp', 'exp','kidsexp','KidsExp'}
        % 
        mymeans=nvec.*exp(problem_all.Invasiveness.kids);
        InvSample=sum(mymeans); 
      case {'InvExp', 'invexp','Invexp', 'invExp'}
        % 
        mymeans=nvec.*exp(0.5*(problem_all.Invasiveness.kids+problem_all.Invasiveness.adults));
        InvSample=sum(mymeans); 
      case {'DRExp', 'drexp','Drexp','DRexp'}
          %
          a= -2;
          b= 0.5;
          myprobs = (1./ (1+exp(-a-b*problem_all.DRscore)));
          newObj = exp(problem_all.Invasiveness.combo).*myprobs; % each genotype's contribution to the total
          InvSample=sum( nvec(:).*newObj(:)); % weighted balance according to how much of each genotype is in the population
  
        
end

InvSample=InvSample/sum(nvec); % carrying capacity. just for scaling. 





