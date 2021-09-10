%% Ellipsoidal membership degree computation

function[w_xk,lambda0,lambda1,f] = mdg(x,m,a,b)
upsilon = zeros(1,size(x,2));
lambda0 = zeros(1,size(x,2));
lambda1 = zeros(1,size(x,2));

for j = 1:size(x,2)
    upsilon(1,j) = ((x(j)-m(j))^2) /  ((b(j)-a(j))^2);
    lambda0(1,j) = (2*(x(j) - m(j))) /  ((b(j)-a(j))^2);
    lambda1(1,j) = (2*(((x(j)-m(j))^2))) / ((b(j)-a(j))^3);
end
theta = 2*sqrt(sum(upsilon));
f = theta / 2;
w_xk = exp(-theta);
end