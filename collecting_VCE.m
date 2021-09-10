Testsets=103;
VCE=zeros(1,Testsets*125000);
Gate_Voltage=VCE;
Collector_current=VCE;
dt=measurement.transient(1).timeDomain.dt  ;
for i=1:Testsets
    if i==1
         VCE(1,(i):125000*i)=measurement.transient(i).timeDomain.collectorEmitterVoltage;
         Gate_Voltage(1,(i):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
          Collector_current(1,(i):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;

    end
   VCE(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.collectorEmitterVoltage;
   Gate_Voltage(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
  Collector_current(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;

end
t=0:dt:(125000*Testsets)*dt-dt;
plot(t,VCE)


for i=1:Testsets
    if i==1
         Gate_Voltage(1,(i):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
         
    end
   Gate_Voltage(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
end

for i=1:Testsets
    if i==1
         Gate_Voltage(1,(i):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
         
    end
   Gate_Voltage(1,(125000*(i-1)+1):125000*i)=measurement.transient(i).timeDomain.gateSignalVoltage;
end