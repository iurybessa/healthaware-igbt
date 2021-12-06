function [rul,rulmax,rulmin,xp]=predictRUL2(eefig,xk,thr,OFFSET)
ngran=numel(eefig);

epsi=1e-3;
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
n=size(xpi,1);
rulvec=inf*ones(size(xpi));
minhi=min(xpi);
% maxhi=max(xpi);
% medianhi=median(xpi);
xd=xpi';
% xd=[xpi,xd(1,1:end-1)];
j=1;
while (minhi<=-epsi && j<=100)
    j=j+1;
    [g,~,~,~] = data_evaluation(eefig,xd,thr);
    if OFFSET
        xold=[1; xd'];
    else
        xold=xd';
    end
    xpi=zeros(size(xk))';
    for i=1:ngran
        xpi=xpi+g(i)*eefig(i).A*xold;
    end
    for k=1:n
        if isinf(rulvec(k)) && xpi(k)<=-epsi
            rulvec(k)=j-1;
        end
    end
    minhi=min(xpi);
    xp(:,j)=xpi;
    xd=[xpi,xd(1,1:end-1)];
end

if j<100
    rul = median(rulvec);
    rulmax = max(rulvec);
    rulmin = min(rulvec);
else
    rul=inf;
    rulmin=inf;
    rulmax=inf;
end

end