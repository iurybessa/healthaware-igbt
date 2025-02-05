function [rul,xp,rul_inf,rul_sup]=predictRUL_v9(eefig,xk,EOL,thr,OFFSET,nu0,rho_nu,zk,pdx,Y)
ngran=numel(eefig);
epsi=0e-3;
if OFFSET
    xd=xk(1,2:end);
else
    xd=xk(1,1:end);
end
[g,~,~,~] = data_evaluation(eefig,zk,thr);
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

xd=[xpi,xd(1,1:end-1)];

j=1;
while ( ((pdx && xp(j,1)<=EOL-epsi) || (~pdx && xp(j,1)>=EOL+epsi) ) && j<=150)
    j=j+1;
    [g,~,~,~] = data_evaluation(eefig,zk,thr);
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
    
    Xi = Hs'*Thetas';
    nuj = Xi*(covmat(nu).*rho_nu)*Xi' + nu0;
    xp_inf(j,1)=xpi-2*vs*sqrt(nuj);
    xp_sup(j,1)=xpi+2*vs*sqrt(nuj);
    
    xd=[xpi,xd(1,1:end-1)];
    nu=[nuj,nu(1:end-1)];
end
if j<150
    rul = j;
    if isempty(find(xp_inf <= EOL, 1)) || find(xp_inf <= EOL, 1)==1
        rul_inf = 0;
        if isempty(find(xp_sup <= EOL, 1))
            rul_sup=0;
        else
            rul_sup = find(xp_sup <= EOL, 1) ;
            rul_inf = rul - (rul_sup-rul);
        end
    else
        rul_inf = find(xp_inf <= EOL, 1) ;
        rul_sup = rul + (rul-rul_inf);
    end
else
    rul=inf;
    rul_inf=inf;
    rul_sup=inf;
end
if rul==rul_inf
    rul;
end
end