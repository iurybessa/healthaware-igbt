%% RUL Estimation with EEFIG

clc, clear all, close all
addpath('EEFIG_FULL');
addpath('to_iury');
%%  Parameters

tau=7; % number of autoregressive terms
OFFSET=0; % If OFFSET=1 then the model has a constant term (bias)
load('features_trig.mat') % IGBT Dataset 
EOL= 1.573; % End of Life
ff=0.999; % forgetting factor
zeta=2; % Required anomalies for creating new rules
buffer=9; % Number of initialization samples (> tau)

%% Initialization

% Data pre-processing
data1 = Mfeatures2(:,2)-EOL;

%data1 = movmean(data1, 10);

idx = hankel(1:tau+1, tau+1:length(data1))';
data2 = data1(idx);
X = data2(:, 1:end-1);
Y = data2(:, end);


idx = hankel(1:tau, tau:length(data1))';
data2 = data1(idx);
X = data2(:, 1:end-1);
Y = data2(:, end);
X = [X (1:length(Y))'];

%%
% RLS initialization
if OFFSET
    Pm0=1e6*eye(tau+1);
    theta{1}=zeros(tau+1,1);
else
	Pm0=1e6*eye(tau);
    theta{1}=zeros(tau,1);
end
P{1}=Pm0;

% EEFIG initialization
[n,p] = size(X);
labels = zeros(n,1);
thr = chi2inv(0.99,p);
separation = 0.2; % c-separation
aux_gran = granule([p,1]);
aux_gran = aux_gran.gran_init(p,X(1:buffer,:));
EEFIG = granule([p,1]);
EEFIG = EEFIG.gran_init(p,X(1:buffer,:));
trackerC = 1*eye(p);
trackerm = mean(X(1:buffer,:));
lambda = 0.99;
Anomalies = [];
continuous_anomalies = 0;

EEFIG_Error = [];
EEFIG_Error_Var = [];


for i=1:buffer
    if OFFSET
        xk = [1 X(i,:)]';
    else
        xk = X(i,:)';
    end
    [thetak,Pk]=wrls_murilo(P{1},Y(i),xk,theta{1},1,1);
    theta{1} = thetak;
    P{1} = Pk;
end

for i = buffer+1:n
    xk = X(i,:);
    %% Change point detection
    [~,~,is_anomaly,~] = data_evaluation(EEFIG,xk,thr);

    %% Change point detection
    [trackerC,trackerm] = tracker_gran(trackerC,trackerm,i,lambda,X(i,:));
    cs = c_separation(EEFIG,trackerC,trackerm,separation);

    %% Anomaly Detection and Data Labeling   
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
 
    %% New EEFIG
    
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
       
    %% Estimation of the A's matrices - Consequent Estimation via RLS
    if OFFSET
        xk = [1 X(i,:)]';
    else
        xk = X(i,:)';
    end
    for k=1:ngran
        [thetak,Pk]=wrls_murilo(P{k},Y(i),xk,theta{k},g(k),1);
        theta{k} = thetak;
        P{k} = Pk;
        EEFIG(k).A=theta{k};
    end
    datahat{i+1} = 0;
    for h = 1:ngran
        datahat{i+1} = datahat{i+1}+g(h)*xk'*theta{h};
    end

    err = Y(i) - datahat{i+1};
    EEFIG_Error = [EEFIG_Error err];

    if i>10
        rho_nu = corr(X)*0+1;
        nu0 = var(EEFIG_Error);
        %nu0 = movvar(EEFIG_Error, 10);
        %nu0 = nu0(end);
        %break
        EEFIG_Error_Var(end+1) = nu0;
        [rul(i-buffer,:),xp]=predictRUL_v4(EEFIG,xk',0,thr,OFFSET,nu0,rho_nu,Y,i);
    else
        rul(i,:) = nan;
    end
end

  

% [g,~,~,~] = data_evaluation(EEFIG,data(i,:),thr)

for i = buffer+2:size(datahat,2)-1
    pred(i,:) = datahat{i};
    deg(i,:) = X(i,1);
end

% c1=colormap(jet(14));
% t=chi2inv(0.999,p);
% for j=1:size(EEFIG,1)
%     hold on;
%     Ellipse_Plot([1/((EEFIG(j).b(1)-EEFIG(j).a(1))/2)^2, 0; 0, 1/((EEFIG(j).b(2)-EEFIG(j).a(2))/2)^2], [EEFIG(j).m],'g');
%     Ellipse_Plot((EEFIG(j).C)/t, EEFIG(j).m,'b');
% end
% grid on
% plot(data(:,1),data(:,2),'k.')

vrul=60-buffer-1:-1:0;
plot(vrul)
hold on
plot(rul)
xlabel('time')
ylabel('RUL')
plot(0.7*vrul,'-.k','Linewidth',2)
plot(1.3*vrul,'-.k','Linewidth',2)
figure
plot(pred)
hold on
plot(deg)


%%
w = X(1:10,:)\Y(1:10,:);
%w = theta{1};
% w=[   -0.2453
%     0.8163
%    -0.9233
%     0.1475
%     0.2123
%    -0.6484
%     1.6484];
% plot(Y)
% hold on
% y = X(1,:)';
% for i = 1:60
%     y(end+1) = y(end-tau+1:end)' * w;
% end
% plot(y(tau+1:end))