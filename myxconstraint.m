function tf = myxconstraint(X,Max)

for n=1:size(X,2)
    thisstring=['mypop(:,' num2str(n) ')= X.s' num2str(n) ';'];
    eval(thisstring);
end

tf = sum(mypop,2)<=Max;
