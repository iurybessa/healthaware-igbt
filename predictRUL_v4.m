function [rul,xp]=predictRUL_v4(eefig,xk,EOL,thr,OFFSET,nu0,rho_nu,Y,k)
ngran=numel(eefig);

if OFFSET
    xd=xk(1,2:end);
else
    xd=xk(1,1:end);
end
[g,~,~,~] = data_evaluation(eefig,xd,thr);
xpi=0;

Thetas = [];
Hs = [];
nu = xd*0;
nu(1)=nu0;

for i=1:ngran
    xpi=xpi+g(i)*xk*eefig(i).A;
    Hs = [Hs; g(i)];
    Thetas = [Thetas eefig(i).A];
end

xp(1,1)=xpi;

vs=1.96;
xp_inf(1,1)=xpi-vs*sqrt(nu0);
xp_sup(1,1)=xpi+vs*sqrt(nu0);

xd=[xd(1,2:end) xpi];
j=1;
while (xp(j,1)<=EOL-1e-3 && j<=100)
    j=j+1;
    [g,~,~,~] = data_evaluation(eefig,xd,thr);
    if OFFSET
        xold=[1 xd];
    else
        xold=xd;
    end
    xpi=0;
    Thetas = [];
    Hs = [];
    for i=1:ngran
        xpi=xpi+g(i)*xold*eefig(i).A;
        Hs = [Hs; g(i)];
        Thetas = [Thetas eefig(i).A];
    end    
    xp(j,1)=xpi;
    
    Xi = Hs'*Thetas((OFFSET+1):end,:)';
    nuj = Xi*(covmat(nu).*rho_nu)*Xi' + nu0;
    xp_inf(j,1)=xpi-vs*sqrt(nuj);
    xp_sup(j,1)=xpi+vs*sqrt(nuj);
    
    %xd=[xpi,xd(1,1:end-1)];
    xd=[xd(1,2:end) xpi];
    nu=[nuj,nu(1:end-1)];
end
1;
if j<100
    rul = j-1;
else
    rul=inf;
end

end