function [pred,deg,rul]=evol_prog(datain,tau,zeta,buffer,ff,sep,EOL,EOLtime)
% This function runs the evolvolving fuzzy modeling algorithm (EEFIG) with
% the "winner takes it all" granularity allocation approach for describe
% the degradation dynamics from a health index (HI) time-series data. The
% input/output variables are described as follows:
%   - datain: column vector with the time-series of the HI;
%   - tau: it is the number of lags in the autoregressive part of the
%   Takagi-Sugeno fuzzy model;
%   - zeta: number of consecutives anomalies necesary to enable the granule
%   creation;
%   - buffer: number of past data used in the RLS;
%   - ff: forgetting factor of the RLS;
%   - sep: c index of the c-saparation metric used for allow the granule
%   creation;
%   - EOL: end of life treshold of the HI to estimate the remaining useful
%   life;
%   - EOLtime: the instant of the EOL.
% 
% Brasilia, October 2021

%% Initialization
data1 = datain';
data = data1(tau:end,1);
for i=1:tau-1
    data=[data,data1(tau-i:end-i,1)];
end
[n,p] = size(data);
theta{1}=zeros(tau+1,1);
Pm0=1e5*eye(tau+1);
P{1}=Pm0;
thr = chi2inv(0.999,p);
separation = sep;

aux_gran = granule([p,1]);
aux_gran = aux_gran.gran_init(p,data(1:buffer,1:p));

EEFIG = granule([p,1]);
EEFIG = EEFIG.gran_init(p,data(1:buffer,1:p));

trackerC = 1*eye(p);
trackerm = mean(data(1:buffer,1:p));
lambda = 0.9;

Anomalies = [];
continuous_anomalies = 0;

%% EEFIG

for i = buffer+1:n
    xk = data(i,:);
    %% Change point detection
    [~,~,is_anomaly,~] = data_evaluation(EEFIG,xk,thr);

    %% Change point detection
    [trackerC,trackerm] = tracker_gran(trackerC,trackerm,i,lambda,data(i,1:p));
    cs = c_separation(EEFIG,trackerC,trackerm,separation);

    %% Anomaly Detection and Data Labeling   
    if (is_anomaly)
        continuous_anomalies = continuous_anomalies + 1;
    else
        continuous_anomalies = 0;
    end

    if is_anomaly > 0
        Anomalies =[Anomalies;xk];
    end
 
    %% New EEFIG
    
    if (cs==1 && continuous_anomalies>(zeta)) && (size(Anomalies,1)-continuous_anomalies+1>0)
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
        if i>buffer
            psi=[];
            for j=1:buffer
                psi_j=[1 data(i-j+1,:)];
                psi=[psi;psi_j];
            end
        else
            psi=[];
            for j=1:i
                psi_j=[1 data(i-j+1,:)];
                psi=[psi;psi_j];
            end
        end
        for k=1:ngran
            theta0=theta{k};
            P0=P{k};
            yk=data(i,1);
            [~,thetap,Pp]=rls_step3(P0,yk,psi,theta0,g(k),ff);
            theta{k}=thetap;
            P{k}=Pp;
            EEFIG(k).A=theta{k};
        end
        datahat{i+1} = 0;
        for h = 1:ngran
            datahat{i+1} = datahat{i+1}+g(h)*[1 data(i,:)]*theta{h};
        end
        if i>10
            [rul(i,:),xp]=predictRUL(EEFIG,[1 data(i,:)],EOL,thr);
        else
            rul(i,:) = nan;
        end     
end

  
%% One-step ahead prediction 
for i = buffer+2:size(datahat,2)-1
    pred(i,:) = datahat{i};
    deg(i,:) = data(i,1);
end
figure(1)
vrul=EOLtime-buffer-1:-1:0;
plot(vrul,'--r','Linewidth',2)
hold on
plot(rul,'b','Linewidth',2)
xlabel('time')
ylabel('RUL')
plot(0.7*vrul,'-.k','Linewidth',2)
plot(1.3*vrul,'-.k','Linewidth',2)
legend('RUL','RUL estimate','Confidence bounds 30%');
title('RUL Prediction');
figure(2)
plot(deg,'r','Linewidth',2)
hold on
plot(pred,'-.k','Linewidth',2)
title('HI regression')
legend('Health index','One-step ahead prediction of HI','Location','southeast')


end