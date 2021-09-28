function [g,basegran,is_anomaly,gran_index] = data_evaluation(basegran,xk,t1)
% Data evaluation
%  Computation of granule membership degree and anomaly condition

sp = 20; % stabilization period
ngran = numel(basegran);
wk = zeros(ngran,1);
lambda0 = zeros(ngran,size(xk,2));
lambda1 = zeros(ngran,size(xk,2));
f = zeros(ngran,1);
for i = 1:ngran
    [wk(i,1),lambda0(i,:),lambda1(i,:),f(i,1)] = mdg (xk,basegran(i).m,...
        basegran(i).a,basegran(i).b);
end
gsum = sum(wk);
g = wk./gsum; % normalizing the membership function value
is_anomaly = 1;
[~,gran_index] = max(g);

for i = 1:ngran
    if(g(i)<1e-6) || (isnan(g(i)))
        continue;
    end
    [basegran(i),alert] = basegran(i).gran_step(xk,sp,t1,wk(i,1),wk,...
        lambda0(i,:),lambda0,lambda1(i,:),lambda1,f(i,1),f,g(i));
    if (alert == 0)
        is_anomaly = 0;
    end
end

if all(~wk)
    gsum = ngran;
    g = ones(ngran,1)./gsum;
else
    gsum = sum(wk);
    g = wk./gsum;
end


% Para avaliar
% for i = 1:ngran
%     [wk(i,1),lambda0(i,:),lambda1(i,:),f(i,1)] = mdg (xk,basegran(i).m,...
%         basegran(i).a,basegran(i).b);
% end
% gsum = sum(wk);
% g = wk./gsum;
% [~,gran_index] = max(g);
end
