%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Nihar Dasari
%Function to calculate droop, settling time and ripple from simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tset Vdroop Vript] =  calctime(load_step_time,Vref);
load('ILoad.mat');
load('Vout.mat');
Tstart = load_step_time;
set = 0;
 ttemp = getsampleusingtime(Vout,Tstart, Tstart+500e-9);
Vdroop = min(ttemp,'Quality',-99,'MissingData','remove');
Vdroop = Vref- Vdroop;
Tset = 1000;
while(set == 0 && Tstart < load_step_time + 500e-9)
    ttemp = getsampleusingtime(Vout,Tstart, Tstart+150e-9);
    Vmax = max(ttemp,'MissingData','remove');
    Vmin = min(ttemp,'MissingData','remove');
    if(Vmax < Vref*1.02 && Vmin > Vref*0.92)
        set = 1;
        Tset = Tstart;
    end
    Tstart = Tstart + 1e-9;
end
Tset = Tset - load_step_time;
 ttemp = getsampleusingtime(Vout,load_step_time+500e-9, load_step_time +750e-9);
 Vript = max(ttemp,'MissingData','remove') - min(ttemp,'MissingData','remove');