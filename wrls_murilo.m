function [w_k, P_k] = wrls_murilo(P_km1, y_k, x_k, w_km1, g, ff)
% Tirei da tese do PAIM

%P_k = P_km1/ff - (g*P_km1*x_k*x_k'*P_km1/ff)/(ff + x_k'*P_km1*x_k);
%w_k = w_km1 + P_k*x_k*g*(y_k - x_k'*w_km1);

a = y_k - x_k'*w_km1;
g = P_km1*x_k/(ff + x_k'*P_km1*x_k);
P_k = (1/ff)*P_km1 - g*x_k'*(1/ff)*P_km1;
w_k = w_km1 + a*g;

end
