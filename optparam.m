% cross-validation

clc, clear all, close all
addpath('EEFIG_FULL');
addpath('data_igbt');
dv=3;
scldfeatname=strcat('device',num2str(dv),'_scaledtrigfeatures.mat')
featname=strcat('device',num2str(dv),'_features.mat')
load(scldfeatname); % IGBT Dataset
load(featname);


tauv=[2:5];
tau2v=[2:5];
ffv=[0.96:0.01:1];
bufferv=[2:6];
zetav=[2:6];

VAL=1; %cross-validation flag
PLOTF=0;
SAVEF=0;    
idcv=0;
bestp=[tauv(1),tau2v(1),ffv(1),bufferv(1),zetav(1)];
bestRA=-1;
for tau=tauv
    for tau2=tau2v
        for ff=ffv
            for buffer=bufferv
                for zeta=zetav
                    idcv=idcv+1
                    clear gvec datahat EEFIG Yk Xk Zk rul vrul pred deg theta Anomalies
                    eefig_prog_v9
                    for kpred=vidx
                        rakp(kpred-vidx(1)+1)=rak(rul(vidx),vrul(vidx)',kpred-vidx(1)+1)*(kpred-vidx(1)+1);
                    end
                    rav(idcv)=sum(rakp);
                    if rav(idcv)>bestRA
                        bestp=[tau,tau2,ff,buffer,zeta];
                        bestRA=rav(idcv);
                    end
                end
            end
        end
    end
end

savename=strcat('bestp_dv',num2str(dv));

save(savename,'bestp')