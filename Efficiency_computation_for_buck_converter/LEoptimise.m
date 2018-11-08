%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Nihar Dasari, PhD student in GREEN laboratory 
%Georgia Institue of Technology.
%Wrapper for efficiency modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [ivropts] = LEoptimise(FSW,L,C,ESR_L,ESR_C,N);
% FSW=91e06; L=18e-9; ESR_L=144e-3; C=22e-9, ESR_C=0.05;
% FSW=119e06; L=12e-9; ESR_L=0.096; C=22e-9, ESR_C=0.05;  
ESR_PCB_PCKG=0.0037; G_PCB_PCKG=0.0178; ESR_Interp=0.0115+5e-3; G_Interp=3.3e-3; off_chip='N'; includeFreqTable='N'; VIN=3.6; tech='28nm'; % off-chip
fstart=1e06;fstep=2e06; fend=FSW;

sizeArray=(fend-fstart)/fstep;
efficiency=zeros(sizeArray,1);

numPhases=1; rTableIndex=2;
Istart=100e-03/numPhases; Imax=1500e-03/numPhases; Istep=20e-03/numPhases;
index3 = 1;
for L = 1e-9:1e-9:40e-9
  %  ESR_L=L*1e9/125;
    for FSW=fstart:fstep:fend
        index2=1+((FSW-fstart)/fstep);
        [e1,l1,pow,ripp]=computeEfficiency(Imax, FSW,L,C,ESR_L,ESR_C, numPhases, rTableIndex,includeFreqTable, VIN, tech, ESR_PCB_PCKG, ESR_Interp, G_Interp, off_chip);
        efficiency(index2,1)=max(e1);
        power = pow(find(max(e1)));
         vr = ripp(find(max(e1)));
%          eff(index2,int16(L/1e-9)) =  efficiency(index2,1);
%          pwr(index2,int16(L/1e-9)) = pow(find(max(e1)));
%          Vrip(index2,int16(L/1e-9)) = vr(find(max(e1)));
        if (vr(find(max(e1))) < 50e-3)
           ivropts(index3).eff =  efficiency(index2,1);  
           ivropts(index3).pow = pow(find(max(e1)));
           ivropts(index3).vr = vr(find(max(e1)));
           ivropts(index3).fsw = FSW;
           ivropts(index3).L = L;
           ivropts(index3).ESR_L = ESR_L;
           index3 = index3 + 1;
        end    
    end
end

  end
            