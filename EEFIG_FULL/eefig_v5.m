%% Evolving Ellipsoidal Fuzzy Information Granules
clc, clear all, close all
%  Initial Parameterization
global r
r = 0;
dec = 10;
tau=5;
%  Data loading
load H1
data1 = L(:,1);
data = data1(tau:end,1);
for i=1:tau-1
    data=[data,data1(tau-i:end-i,1)];
end
Pm0=1e5*eye(tau+1);
P{1}=Pm0;
ff=0.995;
theta{1}=zeros(tau+1,1);
zeta=3;
%  Preprocessing
[n,p] = size(data);
active_gran = zeros(n,1);
labels = zeros(n,1);
thr = chi2inv(0.999,p);
% buffer = p+1;
buffer=2;
LastTen =zeros(10,1);
LastTenIndex = 1;

na = 10;
separation = 0.1;

aux_gran = granule([p,1]);
aux_gran = aux_gran.gran_init(p,data(1:buffer,1:p));

EEFIG = granule([p,1]);
EEFIG = EEFIG.gran_init(p,data(1:buffer,1:p));

trackerC = 1*eye(p);
trackerm = mean(data(1:buffer,1:p));
lambda = 0.9;

active_counter = 0;
active_gran(1:buffer) = ones(buffer,1);
labels(1:buffer) = ones(buffer,1);
Anomalies = [];

continuous_anomalies = 0;
DR = 0;
cond = 0;

w = 1;

A = [];

memb_idx(:,1) = ones(buffer,1);
memb_idx(:,2) = ones(buffer,1);
lb1(:,1) = ones(buffer,1);
lb1(:,2) = ones(buffer,1);
v = 0;
EOL=1.1501;
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
            [K_k,thetap,Pp]=rls_step3(P0,yk,psi,theta0,g(k),ff);
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

  

% [g,~,~,~] = data_evaluation(EEFIG,data(i,:),thr)

for i = buffer+2:size(datahat,2)-1
    pred(i,:) = datahat{i};
    deg(i,:) = data(i,1);
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

