function allvars=makeallvars(numVars)

x1=optimizableVariable('s1',[0 1],'Type','integer');
for n=2:numVars
    mystr=['s' num2str(n)];
    xtemp=optimizableVariable(mystr, [0,1],'Type','integer');
    eval(['x' num2str(n) '=xtemp;'])
end
 xstring = strjoin(cellstr(num2str((1:numVars)', 'x%d')));
 eval(['allvars = [' xstring '];']); 
 