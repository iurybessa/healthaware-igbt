function [rul,xp]=predictRUL(eefig,xk,EOL,thr)
ngran=numel(eefig);

[g,~,~,~] = data_evaluation(eefig,xk,thr);
xpi=zeros(size(xk'));
for i=1:ngran
    xpi=xpi+g(i)*eefig(i).A*xk';
end

xp(1,1)=xpi(end);

j=1;
while (xp(j,1)<=EOL)
    j=j+1;
    [g,~,~,~] = data_evaluation(eefig,xpi',thr);
    xold=xpi;
    xpi=zeros(size(xk'));
    for i=1:ngran
        xpi=xpi+g(i)*eefig(i).A*xold;
    end    
    xp(j,1)=xpi(end);
end

rul = j-1;

end