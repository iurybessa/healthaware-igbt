for j=1:4
    if j==1
        load('Device21.mat')
    elseif j==2
        load('Device31.mat')
    elseif j==3
        load('Device41.mat')
    elseif j==4
        load('Device51.mat')
    end
        
Testsets=90;
VCE=zeros(1,Testsets*125000);
Precursor.Gate_Voltage=VCE;
Precusor.Collector_current=VCE;
dt=measurement.transient(1).timeDomain.dt  ;
% dt=20e-6
% t =(1:length(125000*103))/dt
    features.Mean = zeros(1,Testsets);
    features.Std = zeros(1,Testsets);
    features.Skewness = zeros(1,Testsets);
    features.Kurtosis = zeros(1,Testsets);
    features.Peak2Peak = zeros(1,Testsets);
    features.RMS = zeros(1,Testsets);
    features.CrestFactor = zeros(1,Testsets);
    features.ShapeFactor = zeros(1,Testsets);
    features.ImpulseFactor = zeros(1,Testsets);
    features.MarginFactor = zeros(1,Testsets);
    features.Energy = zeros(1,Testsets);
for i=1:Testsets
    if i==1
         Precursor.VCE(1,(i):125000*i)=measurement.transient(i).timeDomain.collectorEmitterVoltage;
         Precursor.Gate_Voltage(1,(i):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
         Precursor.Collector_current(1,(i):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
    features.Mean(i) = mean(VCE(1,(i):125000*i));
    features.Std(i) = std(VCE(1,(i):125000*i));
    features.Skewness(i) = skewness(VCE(1,(i):125000*i));
    features.Kurtosis(i) = kurtosis(VCE(1,(i):125000*i));
    features.Peak2Peak(i) = peak2peak(VCE(1,(i):125000*i));
    features.RMS(i) = rms(VCE(1,(i):125000*i));
    features.CrestFactor(i) = max(VCE(1,(i):125000*i))/features.RMS(i);
    features.ShapeFactor(i) = features.RMS(i)/mean(abs(VCE(1,(i):125000*i)));
    features.ImpulseFactor(i) = max(VCE(1,(i):125000*i))/mean(abs(VCE(1,(i):125000*i)));
    features.MarginFactor(i) = max(VCE(1,(i):125000*i))/mean(abs(VCE(1,(i):125000*i)))^2;
    features.Energy(i)= sum(VCE(1,(i):125000*i).^2);
    end
   VCE(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.collectorEmitterVoltage;
   Gate_Voltage(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
  Collector_current(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
    features.Mean(i) = mean(VCE(1,(125000*(i-1)+1):125000*i));
    features.Std(i) = std(VCE(1,(125000*(i-1)+1):125000*i));
    features.Skewness(i) = skewness(VCE(1,(125000*(i-1)+1):125000*i));
    features.Kurtosis(i) = kurtosis(VCE(1,(125000*(i-1)+1):125000*i));
    features.Peak2Peak(i) = peak2peak(VCE(1,(125000*(i-1)+1):125000*i));
    features.RMS(i) = rms(VCE(1,(125000*(i-1)+1):125000*i));
    features.CrestFactor(i) = max(VCE(1,(125000*...
        (i-1)+1):125000*i))/features.RMS(i);
    features.ShapeFactor(i) = features.RMS(i)/mean(abs...
        (VCE(1,(125000*(i-1)+1):125000*i)));
    features.ImpulseFactor(i) = max(VCE(1,(125000*(i-1)+1):125000*i))...
        /mean(abs(VCE(1,(125000*(i-1)+1):125000*i)));
    features.MarginFactor(i) = max(VCE(1,(125000*...
        (i-1)+1):125000*i))/mean(abs(VCE(1,(125000*(i-1)+1):125000*i)))^2;
    features.Energy(i) = sum(VCE(1,(125000*(i-1)+1):125000*i).^2);
end
t=0:dt:(125000*Testsets)*dt-dt;
figure(1)
plot(t,VCE) 
figure(2)
subplot(2,2,1)
plot(features.Energy)
title('Energy')
subplot(2,2,2)
plot(features.Mean)
title('mean')
subplot(2,2,3)
plot(features.Std);
title('std')
subplot(2,2,4)
plot(features.Skewness);
title('skewness')
figure(3)
subplot(2,2,1)
plot(features.Kurtosis);
title('kurtosis')
subplot(2,2,2)
plot(features.Peak2Peak);
title('Peak2peak')
subplot(2,2,3)
plot(features.RMS);
title('RMS')
subplot(2,2,4)
plot(features.CrestFactor);
title('Crestfactor')
figure(4)
subplot(3,1,1)
plot(features.ShapeFactor);
title('shapefactor')
subplot(3,1,2)
plot(features.ImpulseFactor);
title('shapefctor')
subplot(3,1,3)
plot(features.MarginFactor);
title('marginfactor')
W_RMS=sphx(features.RMS,'PCA');
hold on
plot(W_RMS)
plot(features.RMS)

%% *********************************
% [U,S,V]= svd(features.CrestFactor','econ');
% % Least-square fit
% xtilde = V*inv(S)*U'*features.CrestFactor';
% plot(xtilde*features.CrestFactor','b--')
X=[features.RMS',features.Kurtosis',features.Mean',features.Skewness',features.Peak2Peak',...
    features.CrestFactor',features.Std',features.Energy'];
M=mean(X);
STD=std(X);
Normalized = (X - M)./STD;
[L,A]=RPCA(Normalized);
coef=pca(Normalized);
PCA1 = (X - M) ./ STD * coef(:, 1);
PCA2 = (X - M) ./ STD * coef(:, 2);

figure(5)
hold on
if j==1
plot(L(:,1),'b-*')

elseif j==2
plot(L(:,1),'k-*')
elseif j==3
plot(L(:,1),'r-o')
elseif j==4
plot(L(:,1),'k-s')

legend({'device1','device2','device3','device4'},'location','northwest')
xlabel('cycles')
ylabel('healthindex')

end
end