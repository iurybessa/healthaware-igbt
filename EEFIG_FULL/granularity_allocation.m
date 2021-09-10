%% Granularity allocation

function [updated_gran,is_improved]=granularity_allocation(gran,gran_test,mahal_xk,g,wk_i,wk,l0_i,l0,l1_i,l1,f_i,f)
p = size(l1,2);

aux_l0 = zeros(numel(wk),p);
aux_l1 = zeros(numel(wk),p);

sum_w = sum(wk);
sum_w2 = (sum_w)^2;

sum_l0 = zeros(1,p);
sum_l1 = zeros(1,p);

for i=1:numel(wk)
    for j = 1:p
        aux_l0(i,j) = (wk(i,1) * l0(i,j)) / (f(i,1));
        aux_l1(i,j) = (wk(i,1) * l1(i,j)) / (f(i,1));
    end
end

for j = 1:p
    sum_l0(1,j) = sum(aux_l0(:,j));
    sum_l1(1,j) = sum(aux_l1(:,j));
end

dg_da = ((wk_i * sum_l1) / (sum_w2)) - ((wk_i * l1_i) / (f_i * sum_w));
dg_dm = ((wk_i * l0_i) / (f_i * sum_w)) - ((wk_i * sum_l0) / (sum_w2));
dg_db = ((wk_i * l1_i) / (f_i * sum_w)) - ((wk_i * sum_l1) / (sum_w2));

delta_a = gran.a - gran_test.a;
delta_m = gran.m - gran_test.m;
delta_b = gran.b - gran_test.b;

gran.gsum = gran.gsum + (g - (dg_da * delta_a') + (dg_dm * delta_m') ...
    + (dg_db * delta_b'));

Q_it = (mahal_xk * gran.gsum);
gran.Q_it = [gran.Q_it;Q_it];
gran.Q = cumsum(gran.Q_it);

if (gran.Q(end)) >  (gran_test.Q(end))
    updated_gran = gran;
    is_improved = 1;
    updated_gran.is_improved = is_improved;
    updated_gran.counter = updated_gran.counter + 1;
else
    updated_gran = gran_test;
    is_improved = 0;
    updated_gran.is_improved = is_improved;
end
end