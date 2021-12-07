clc; clear all; close all;
addpath('EEFIG_FULL');
addpath('to_iury');

%%  Parameters
tau = 3; % number of autoregressive terms
ff = 0.999; % forgetting factor
zeta = 2; % Required anomalies for creating new rules
buffer = 4; % Number of initialization samples (> tau)
separation = 0.2; % c-separation
lambda = 0.95;

%% Initialization
load('features_trig.mat') % IGBT Dataset 
EOL= 1.573; % End of Life

% Data pre-processing
data1 = Mfeatures2(:,2)-EOL;

%data1 = movmean(data1, 10);

idx = hankel(1:tau+1, tau+1:length(data1))';
data2 = data1(idx);
X = data2(:, 1:end-1);
Y = data2(:, end);

[n,p] = size(X);
thr = chi2inv(0.99,p);

%%
% RLS initialization
theta{1}=zeros(p,1);
P{1}=1e6*eye(p);

% EEFIG initialization
labels = zeros(n,1);
aux_gran = granule([p,1]);
aux_gran = aux_gran.gran_init(p,X(1:buffer,:));
EEFIG = granule([p,1]);
EEFIG = EEFIG.gran_init(p,X(1:buffer,:));
trackerC = 1*eye(p);
trackerm = mean(X(1:buffer,:));
trackerC = inv(cov(X(1:buffer,:)));
EEFIG.C = trackerC;
Anomalies = [];
continuous_anomalies = 0;

EEFIG_Error = [];
EEFIG_Error_Var = [];

%%
for i = buffer+1:n
    xk = X(i,:);
    % Change point detection
    [~,~,is_anomaly,~] = data_evaluation(EEFIG,xk,thr);

    % Change point detection
    [trackerC,trackerm] = tracker_gran(trackerC,trackerm,i,lambda,X(i,:));
    cs = c_separation(EEFIG,trackerC,trackerm,separation);

    % Anomaly Detection and Data Labeling   
    if (is_anomaly)
        continuous_anomalies = continuous_anomalies + 1;
    else
        continuous_anomalies = 0;
    end

    if is_anomaly > 0
        Anomalies =[Anomalies;xk];
    else
        Anomalies=[];
    end
 
    % New EEFIG
    
    if (cs==1 && continuous_anomalies>(zeta))
        newEEFIG = aux_gran.gran_init(p,Anomalies);
        Anomalies = [];
        EEFIG = [EEFIG;newEEFIG];
        ngran = numel(EEFIG);
        P{ngran}=Pm0;
        theta{ngran}=theta{ngran-1};
    end

    [g,EEFIG,~,lastactive] = data_evaluation(EEFIG,xk,thr);
    ngran = numel(EEFIG);   
       
    clf;
    hold on
    plot3(X(1:i,1), X(1:i,2), X(1:i,3), 'ko')
    Ellipse_Plot(EEFIG.C/chi2inv(0.99,p), EEFIG.m, 'r')
    Ellipse_Plot(trackerC/chi2inv(0.99,p), trackerm, 'b')
    hold off
    pause(0.01);
end

