function [rul,xp]=predictRUL2(eefig,xk,EOL,thr,OFFSET)
ngran=numel(eefig);

if OFFSET
    xd=xk(1,2:end);
else
    xd=xk(1,1:end);
end
[g,~,~,~] = data_evaluation(eefig,xd,thr);
xpi=zeros(size(xk))';
for i=1:ngran
    xpi=xpi+g(i)*eefig(i).A*xk';
end

xp(:,1)=xpi;
xd=[xpi,xd(1,1:end-1)];
j=1;
while (xp(j,1)<=EOL && j<=100)
    j=j+1;
    [g,~,~,~] = data_evaluation(eefig,xd,thr);
    if OFFSET
        xold=[1 xd];
    else
        xold=xd;
    end
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