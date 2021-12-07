clc; clear all; close all;

%%  Parameters
ff = 1; % forgetting factor
tau = 3;
tp = 5;

%% Initialization
load('features_trig.mat') % IGBT Dataset 
EOL = 1.573; % End of Life
data = Mfeatures2(:,2) - EOL;
%data = log(data-min(data)+eps);
%data = data(2:end);
%data = data-max(data);
%data = smoothdata(data, 'lowess');
%data = movmean(data, 10);

idx = hankel(1:tau+1, tau+1:length(data))';
data = data(idx);
X = data(:, 1:end-1);
Y = data(:, end);

%%
[n, D] = size(X);

w = zeros(D,1);
P = eye(D)*1e5;

for i = 1:60
    xk = X(i,:)';
    yk = Y(i);
    
    % Estimation
    y(i) = xk'*w;
    
    YY = [X(1,:)'; Y];
    dy(i) = mean(diff(YY(max(1,i-5):i)));
    if dy(i) > 0.0025
    
    % Parameter update RLS
    a = yk - y(i);
    g = P*xk/(ff + xk'*P*xk);
    P = (1/ff)*P - g*xk'*(1/ff)*P;
    w = w + a*g;
    
    end
    
    % Prognostics
    if i >= tp
        f = Y(1:i);
        
        for k = 1:length(Y)-i
            xk = f(end-tau+1:end);
            f(end+1) = xk'*w;
        end
    
        f(1:i-1) = NaN;
        
        idx = find(f > -1e-2, 1);
        
        clf;
        hold on
        plot(Y);
        plot(f);
        plot(i,Y(i), 'ko');
        plot(idx,Y(idx), 'k*');
        hold off;
        
        if isempty(idx)
            rul(i) = Inf;
        else
            rul(i) = idx - i;
        end
        
        pause(0.1);
    end
end


%%
figure;
vrul = 60-tau:-1:0;
hold on
plot(vrul, 'b-')
plot(vrul*1.3, 'b:')
plot(vrul*0.7, 'b:')
plot(rul, 'r--')
hold off
