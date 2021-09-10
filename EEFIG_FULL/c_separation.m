%% CSeparation condition function

function [c_separated] = c_separation(basegran,c2,m2,c)
c_separated = 1;
if (sum(sum(isnan(c2)))>0)
    return;
end

if (sum(sum(isinf(c2)))>0)
    return;
end

for j=1:size(basegran,1)
    temp1 = norm(basegran(j).m-m2);
    T1 = eig(basegran(j).C);
    T2 = eig(c2);
    temp2 = c*sqrt(numel(c2)*max([1/((T1(1))) 1/((T2(1)))]));
%     temp2 = c*sqrt(numel(m2)*max([max(eig((basegran(j).C))),max(eig((c2)))]));

    if (temp1<temp2)
        c_separated = 0;
    end
end
end