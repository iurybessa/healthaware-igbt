clc; clear all; close all;

%%  Parameters
ff = 0.999; % forgetting factor
tau = 1;
tp = 30;

%% Initialization
load('features_trig.mat') % IGBT Dataset 
EOL = 1.573; % End of Life
data = Mfeatures2(:,2) - EOL;
%data = movmean(data, 10);

idx = hankel(1:tau+1, tau+1:length(data))';
data = data(idx);
X = data(:, 1:end-1);
X = [sqrt((1:size(X,1))'/1000) X];
Y = data(:, end);

%%
[n, D] = size(X);

w = zeros(D,1);
P = eye(D)*1e5;

for i = 1:n
    xk = X(i,:)';
    yk = Y(i);
    
    % Estimation
    y(i) = xk'*w;
    
    a = yk - y(i);
    g = P*xk/(ff + xk'*P*xk);
    P = (1/ff)*P - g*xk'*(1/ff)*P;
    w = w + a*g;
    
    f = Y(1:i);
    if i >= tp
        for k = 1:80
            xk = [sqrt((i+k)/1000);f(end-tau+1:end)];
            f(end+1) = xk'*w;
        end

        hold on
        plot(Y);
        plot(f);
        hold off;
        break
    end
end