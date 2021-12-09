%% RUL Estimation with EEFIG

clc, clear all, close all
addpath('EEFIG_FULL');
addpath('data_igbt');
%%  Parameters

tau=10; % number of autoregressive terms
tau2=9; % number of autoregressive terms
OFFSET=0; % If OFFSET=1 then the model has a constant term (bias)
load('device2_scaledtrigfeatures.mat') % IGBT Dataset
load('device2_features.mat')
iEOL=68;
iFeat=4;
pdx=0;
EOL= Mfeatures2(iEOL,iFeat); % End of Life
ff=0.999; % forgetting factor
zeta=3; % Required anomalies for creating new rules
buffer=5; % Number of initialization samples (> tau)

%% Initialization

% Data pre-processing
%EOL = Mfeatures2(60,4);
data1 = Mfeatures2(:,iFeat)-EOL;
idx = hankel(1:tau+1, tau+1:length(data1))';
idx2 = hankel(1:tau2, tau2:length(features.Energy))';
data2 = data1(idx);
Xk = data2(:, 1:end-1);
Yk = data2(:, end);
% Xk = [Xk (1:length(Yk))'];
Zk = features.Energy(idx2);


% RLS initialization
if OFFSET
    Pm0=0.1e3*eye(tau+1);
    theta{1}=zeros(tau+1,1);
else
	Pm0=0.1e3*eye(tau);
    theta{1}=zeros(tau,1);
end
P{1}=Pm0;
for j=1:buffer
    theta0=theta{1};
    P1=P{1};
    [~,t2,P2]=rls_step3(P1,Yk(j),Xk(j,:),theta0,1,1);
%             [t2,P2]=wrls_murilo(P1,Yk(j),Xk(j,:)',theta0,g(k),ff);
    P{1}=P2;
    theta{1}=t2;
end

% EEFIG initialization
[n,p] = size(Zk);
labels = zeros(n,1);
thr = chi2inv(0.999,p);
separation = 2; % c-separation
aux_gran = granule([p,1]);
aux_gran = aux_gran.gran_init(p,Zk(1:buffer,1:p));
EEFIG = granule([p,1]);
EEFIG = EEFIG.gran_init(p,Zk(1:buffer,1:p));
trackerC = 1*eye(p);
trackerm = mean(Zk(1:buffer,1:p));
lambda = 0.9;
Anomalies = [];
continuous_anomalies = 0;

EEFIG_Error = [];
EEFIG_Error_Var = [];

for i = buffer+1:n-tau
    xk = Zk(i,:);
    %% Change point detection
    [~,~,is_anomaly,~] = data_evaluation(EEFIG,xk,thr);

    %% Change point detection
    [trackerC,trackerm] = tracker_gran(trackerC,trackerm,i,lambda,Zk(i,:));
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
        continuous_anomalies=0;
        EEFIG = [EEFIG;newEEFIG];
        ngran = numel(EEFIG);
%         P{ngran}=Pm0;
        P{ngran}=P{ngran-1};
        theta{ngran}=theta{ngran-1};
%         theta{ngran}=(Yk(i-zeta+1:i)\Xk(i-zeta+1:i,:))';
        for j=i-zeta:i
            theta0=theta{ngran};
            P1=P{ngran};
            [~,t2,P2]=rls_step3(P1,Yk(j),Xk(j,:),theta0,g(k),1);
%             [t2,P2]=wrls_murilo(P1,Yk(j),Xk(j,:)',theta0,g(k),ff);
            P{ngran}=P2;
            theta{ngran}=t2;
        end
%         P{ngran}=inv(cov(Xk(i-zeta+1:i,:)));
%         theta{ngran}=(Yk(i-zeta+1:i)\Xk(i-zeta+1:i,:))';
%         theta{ngran}=t2;
%         theta{ngran-1};
        gvec=[gvec,zeros(size(gvec,1),1)];       
    end

    [g,EEFIG,~,lastactive] = data_evaluation(EEFIG,xk,thr);
    gvec(i,:)=g;
    ngran = numel(EEFIG);   
       
    %% Estimation of the A's matrices - Consequent Estimation via RLS
%         if i>buffer
%             psi=[];
%             for j=1:buffer
%                 if OFFSET
%                     psi_j=[1 data(i-j,:)];
%                 else
%                     psi_j=[data(i-j,:)];
%                 end
%                 psi=[psi;psi_j];
%             end
%         else
%             psi=[];
%             for j=1:i
%                 if OFFSET
%                     psi_j=[1 data(i-j,:)];
%                 else
%                     psi_j=[data(i-j,:)];
%                 end
%                 psi=[psi;psi_j];
%             end
%         end        
        for k=1:ngran
            theta0=theta{k};
            P0=P{k};
%             yk=flipud(data(i-buffer+1:i,1));
            [K_k,thetap,Pp]=rls_step3(P0,Yk(i),Xk(i,:),theta0,g(k),ff);
            theta{k}=thetap;
            P{k}=Pp;
            EEFIG(k).A=theta{k};
        end
        datahat{i+1} = 0;
        ykh = 0;
        for h = 1:ngran
            if OFFSET
                datahat{i+1} = datahat{i+1}+g(h)*[1 data(i,:)]*theta{h};
            else
                datahat{i+1} = datahat{i+1}+g(h)*[Xk(i,:)]*theta{h};
                ykh = ykh+g(h)*[Xk(i,:)]*theta{h};
            end
        end
        
        err = Yk(i) - ykh;
        EEFIG_Error = [EEFIG_Error err];
        
        if i>20
            if OFFSET
                [rul(i,:),xp]=predictRUL(EEFIG,[1 data(i,:)],EOL,thr,OFFSET);
            else
                %triu([0.999.^(0:2); 0.999.^(-1:1); 0.999.^(-2:0)])
                rho_nu = corr(Xk)*0+1;
                
                nu0 = var(EEFIG_Error(max(1,end-buffer+1):end));
                EEFIG_Error_Var(end+1) = nu0;
%                 if Yk(i)<=0
%                     i
%                 end
                [rul(i-buffer,:),xp,rul_inf(i-buffer,:),rul_sup(i-buffer,:)]=predictRUL(EEFIG,[Xk(i,:)],0,thr,OFFSET,nu0,rho_nu,xk,pdx,Yk(i:end));
            end
        else
            rul(i,:) = nan;
        end
end

  

% [g,~,~,~] = data_evaluation(EEFIG,data(i,:),thr)

for i = buffer+2:size(datahat,2)-1
    pred(i,:) = datahat{i};
    deg(i,:) = Xk(i,1);
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

%%
subplot(4,1,1:3)
vrul=iEOL-buffer-tau-1:-1:0;

vidx = 1:length(vrul);
plot(vidx,vrul(vidx),'b--','Linewidth',1)
hold on
plot(vidx,rul(vidx),'r-','LineWidth',1.5)
%plot(vidx,rul_inf(vidx))
% rul_sup = rul + (rul-rul_inf);
%plot(vidx,rul_sup(vidx));
errorbar(vidx,rul(vidx),rul(vidx)-rul_inf(vidx),rul(vidx)-rul_sup(vidx),'ko')

xlabel('time')
ylabel('RUL')
alpha=0.3;
plot(vidx,(1-alpha)*vrul(vidx),':b','Linewidth',1.5)
plot(vidx,(1+alpha)*vrul(vidx),':b','Linewidth',1.5)
[~,maxgvec] = max(gvec');

subplot(4,1,4)
stairs(vidx, maxgvec(vidx)*5, 'LineWidth', 2);


%%
figure
plot(pred)
hold on
plot(deg)
figure
plot(gvec)
legend('g_1','g_2','g_3','g_4','g_5','g_6')

