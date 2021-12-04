function [m] = covmat(vs)
%COVMAT Summary of this function goes here
%   Detailed explanation goes here
S = length(vs);
m = zeros(S);
for i=1:S
    m(i,i) = vs(i);
    for j = i+1:S
        m(i,j) = sqrt(vs(i))*sqrt(vs(j));
    end
end
m = m+tril(m',-1);

end

