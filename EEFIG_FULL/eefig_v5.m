%% Evolving Ellipsoidal Fuzzy Information Granules
clc, clear all, close all
%  Initial Parameterization
global r
r = 0;
dec = 10;

%  Data loading
load ruspini
data1 = ruspini;
data = data1(:,1:2);

%  Preprocessing
[n,p] = size(data);
active_gran = zeros(n,1);
labels = zeros(n,1);
thr = chi2inv(0.999,p);
buffer = p+1;
LastTen =zeros(10,1);
LastTenIndex = 1;

na = 10;
separation = 2;

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
Qs=zeros(buffer,1);
for i = buffer+1:n
    %% Primary Level
    [g,EEFIG,is_anomaly,lastactive] = data_evaluation(EEFIG,data(i,:),thr);
    for i1 = 1 : numel(EEFIG)
        Q{i1}(i) = EEFIG(i1).Q(end);
%         imp{i1}(i) = EEFIG(i1).is_improved;
    end
    %% Change point detection
    [trackerC,trackerm] = tracker_gran(trackerC,trackerm,i,lambda,data(i,1:p));
    cs = c_separation(EEFIG,trackerC,trackerm,separation);
    LastTen(LastTenIndex) = is_anomaly;
    LastTenIndex = LastTenIndex + 1;
    %% Anomaly Detection and Data Labeling
    if (is_anomaly)
        continuous_anomalies = continuous_anomalies + 1;
    else
        continuous_anomalies = 0;
    end
    if LastTenIndex > 10
        LastTenIndex = 1;
    end
    labels(i)=lastactive;
    if (lastactive == active_gran(i-1,1))
        active_gran(i,1) = lastactive;
        active_counter = na;
    else
        if (active_counter > 0)
            active_gran(i,1) = active_gran(i-1,1);
        else
            active_gran(i,1) = lastactive;
            active_counter = na+1;
        end
        active_counter = active_counter-1;
    end
    
    if is_anomaly > 0
        Anomalies =[Anomalies;data(i,:)];
    end
    %% New EEFIG
    if (cs==1 && continuous_anomalies>(p)) && (size(Anomalies,1)-continuous_anomalies+1>0)
        %         if (cs==1 && continuous_anomalies>(p))
%         cond(i,1) = 1;
        newEEFIG = aux_gran.gran_init(p,Anomalies(size(Anomalies,1)...
            -continuous_anomalies+1:size(Anomalies,1),1:p));
        labels(i-continuous_anomalies+1:i)=ones(continuous_anomalies,1).*numel(EEFIG);
        Anomalies(size(Anomalies,1)-continuous_anomalies+1:size(Anomalies,1),:) = [];
        %                 Anomalies = [];
        EEFIG = [EEFIG;newEEFIG];
        active_gran(i-9:i,1) = ones(10,1)*numel(EEFIG);
    end
    %% Updating the membership degrees with the new base of granules
    ngran = size(EEFIG,1);
    [g,~,~,~] = data_evaluation(EEFIG,data(i,:),thr);
    act{i,1} = g;
    [memb_idx(i,1),memb_idx(i,2)] = max(act{i,1});
    
    %% Estimation of the A's matrices - Consequent E////stimation via RLS
    if i>21
        w = 20;
        psik = [];
        Ak = [];
        for j = 1:p
            xr{j} = data(i-w:i,j);
            psik = [psik,data(i-w-1:i-1,j)];
        end
        for j = 1:p
            theta{j} = pinv(psik)*xr{j};
            Ak = [Ak;theta{j}'];
        end
        for l=1:ngran
            A{i,l}=[];
            if l==memb_idx(i,2)
                A{i,l}=Ak;
            elseif isempty(A{i-1,l})
                A{i,l}=Ak;
            else
                A{i,l}=A{i-1,l};
            end
        end
        datahat{i+1} = 0;
        for h = 1:ngran
            datahat{i+1} = datahat{i+1}+act{i}(h)*A{i,h}*data(i,:)';
        end
        
    end
    
    for j=1:ngran
        Qs(i,j)=log(EEFIG(j).Q(end)+1);
    end
end


for i = 23:size(datahat,2)
    pred(i,:) = datahat{i};
end

c1=colormap(jet(14));
t=chi2inv(0.999,p);
for j=1:size(EEFIG,1)
    hold on;
    Ellipse_Plot([1/((EEFIG(j).b(1)-EEFIG(j).a(1))/2)^2, 0; 0, 1/((EEFIG(j).b(2)-EEFIG(j).a(2))/2)^2], [EEFIG(j).m],'g');
    Ellipse_Plot((EEFIG(j).C)/t, EEFIG(j).m,'b');
end
grid on
plot(data(:,1),data(:,2),'k.')

