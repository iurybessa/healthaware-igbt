function [rul,xp]=predictRUL(eefig,xk,EOL,thr)
ngran=numel(eefig);

xd=xk(1,2:end);
[g,~,~,~] = data_evaluation(eefig,xd,thr);
xpi=0;
for i=1:ngran
    xpi=xpi+g(i)*xk*eefig(i).A;
end

xp(1,1)=xpi;
xd=[xpi,xd(1,1:end-1)];
j=1;
while (xp(j,1)<=EOL && j<=100)
    j=j+1;
    [g,~,~,~] = data_evaluation(eefig,xd,thr);
    xold=[1 xd];
    xpi=0;
    for i=1:ngran
        xpi=xpi+g(i)*xold*eefig(i).A;
    end    
    xp(j,1)=xpi;
    xd=[xpi,xd(1,1:end-1)];
end

if j<100
    rul = j-1;
else
    rul=inf;
end

end