%% RUL Estimation with EEFIG

clc, clear all, close all
addpath('EEFIG_FULL');
addpath('to_iury');
%%  Parameters

tau=4; % number of autoregressive terms
OFFSET=0; % If OFFSET=1 then the model has a constant term (bias)
load('features_trig.mat') % IGBT Dataset 
EOL=1.573; % End of Life
ff=0.999; % forgetting factor
zeta=2; % Required anomalies for creating new rules
buffer=5; % Number of initialization samples (> tau)

%% Initialization

% Data pre-processing
datain=Mfeatures2(:,2:end);
datain(:,3)=1./datain(:,3);
datain(:,6)=1./datain(:,6);
datain(:,8)=1./datain(:,8);
datain(:,9)=1./datain(:,9);
datain(:,10)=1./datain(:,10);
% datain(:,11)=1./datain(:,11);
scale=kron(datain(60,:),ones(size(datain,1),1));
D1=datain./scale -1;
data1=[D1(:,1:3),D1(:,6:8),D1(:,11)];
nfeat=size(data1,2);
data = data1(tau:end,:);
for i=1:tau-1
    data=[data,data1(tau-i:end-i,:)];
end

% RLS initialization
for l=1:nfeat
    if OFFSET
        Pm0=1e5*eye(nfeat*(tau+1));
        theta{1,l}=zeros(nfeat*(tau+1),1);
    else
        Pm0=1e5*eye(nfeat*tau);
        theta{1,l}=zeros(nfeat*tau,1);
    end
    P{1,l}=Pm0;
end

% EEFIG initialization
[n,p] = size(data);
labels = zeros(n,1);
thr = chi2inv(0.99,p);
separation = 0.2; % c-separation
aux_gran = granule([p,1]);
aux_gran = aux_gran.gran_init(p,data(1:buffer,1:p));
EEFIG = granule([p,1]);
EEFIG = EEFIG.gran_init(p,data(1:buffer,1:p));
trackerC = 1*eye(p);
trackerm = mean(data(1:buffer,1:p));
lambda = 0.9;
Anomalies = [];
continuous_anomalies = 0;

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
    else
        Anomalies=[];
    end
 
    %% New EEFIG
    
    if (cs==1 && continuous_anomalies>(zeta))
        newEEFIG = aux_gran.gran_init(p,Anomalies);
        Anomalies = [];
        EEFIG = [EEFIG;newEEFIG];
        ngran = numel(EEFIG);
        for l=1:nfeat
            P{ngran,l}=Pm0;
            theta{ngran,l}=theta{ngran-1,l};
        end
    end

    [g,EEFIG,~,lastactive] = data_evaluation(EEFIG,xk,thr);
    ngran = numel(EEFIG);   
       
    %% Estimation of the A's matrices - Consequent Estimation via RLS
        if i>buffer
            psi=[];
            YK=[];
            for j=1:buffer
                if OFFSET
                    psi_j=[1 data(i-j,:)];
                else
                    psi_j=[data(i-j,:)];
                end
                YK=[YK;data(i-j+1,1:nfeat)];
                psi=[psi;psi_j];
            end
        else
            psi=[];
            for j=1:i
                if OFFSET
                    psi_j=[1 data(i-j,:)];
                else
                    psi_j=[data(i-j,:)];
                end
                YK=[YK;data(i-j+1,1:nfeat)];
                psi=[psi;psi_j];
            end
        end
        for k=1:ngran
            if OFFSET
                Ak=zeros(nfeat,nfeat*(tau+1));
            else
                Ak=zeros(nfeat,nfeat*tau);
            end
            for l=1:nfeat
                P0=P{k,l};
                theta0=theta{k,l};
                yk=YK(:,l);
                [K_k,thetap,Pp]=rls_step3(P0,yk,psi,theta0,g(k),ff);
                theta{k,l}=thetap;
                Ak(l,:)=thetap;
                P{k,l}=Pp;
            end            
%             theta0=theta{k};
%             P0=P{k};
            EEFIG(k).A=Ak;
        end
        datahat{i+1} = 0;
        for h = 1:ngran
            if OFFSET
                datahat{i+1} = datahat{i+1}+g(h)*EEFIG(h).A*psi(1,:)';
            else
                datahat{i+1} = datahat{i+1}+g(h)*EEFIG(h).A*psi(1,:)';
            end
        end
%         if i>10
%             if OFFSET
%                 [rul(i,:),rulmax(i,:),rulmin(i,:),xp]=predictRUL2(EEFIG,psi(1,:),thr,OFFSET);
%             else
%                 [rul(i,:),rulmax(i,:),rulmin(i,:),xp]=predictRUL2(EEFIG,psi(1,:),thr,OFFSET);
%             end
%         else
%             rul(i,:) = nan;
%             rulmin(i,:) = nan;
%             rulman(i,:) = nan;
%         end     
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

