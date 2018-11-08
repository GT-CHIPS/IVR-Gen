%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Nihar Dasari, PhD student in GREEN laboratory 
%Georgia Institue of Technology.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Usage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the target specs from target_specs.csv
% Format of the target_specs.csv follows
%
% Vin                   Input voltage in volts
% V_ref                 Reference voltage in volts
% ref_step_en           (1/0) Enables/Disables reference transient event
% ref_step              Reference step in V to reach Vref (Used for reference transient simulations)
% F_SW                  Power stage/DPWM switching Frequency in Hz
% N                     Factor of F_SW by which ADC samples the voltage error. N!=1 for multisampling
% phm_d                 Target phase margin in degrees
% Fc                    Target Unity gain or crossover frequency in Hz
% L                     Output filter inductance
% C                     Output filter capacitance
% ESR_L                 Effective series resistance of output filter inductor
% ESR_C                 Effective series resistance of output filter capacitor
% L_BW                  Inducatance for the pads
% ESR_L_BW              Effective series resistance for PAD inductance
% C_DECAP               Input decap at pads
% ESR_C_DECAP           Effecrtive series resistance for input decaps at pads
% I_Load                Load current
% ADC_lower_range       Lower range of ADC
% ADC_higher_range      +Z_ovrUpper range of ADC
% ADC_reso_bits         ADC resolution in bits
% DPWM_reso_bits        DPWM resolution in bits
% load_step_en          (1/0) Enables/Disables load transient event
% I_load_init           Initial output current in A before load transient event (Used for load transient simulations)
% load_step_time        Time in s at which load step is introduced (Used for load transient simulations)
% load_step             Load step value in A at load_step_time (Used for load transient simulations)
% wk0_bits
% fb_wk1_bits
% fb_wk2_bits
% fb_wk3_bits
% comp_gain             Compensator gain
% quiet                 (1/0) Disables/Enables detailed log statements
% simu                  Simulink Simulation
% auto                  Automatically optimise for Bandwidth and phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
% close all
addpath('./Efficiency_computation_for_buck_converter')
% IVR
%%% Importing target specs from csv file %%%%%
spec_file_name='target_specs_ck.csv';
fileID=fopen(spec_file_name,'r');
target_specs= textscan(fileID,'%s%f','Delimiter',',');
fclose(fileID);
spec_names=matlab.lang.makeValidName(target_specs{1});
spec_vals=target_specs{2};

for i=1:length(spec_names)
    evalc([spec_names{i} ' = spec_vals(i)']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Transient event conditions %%%%%%%%%%%
if (ref_step_en==0)
    ref_step=0;
end
vref=V_ref-ref_step;


if (load_step_en==0)
    load_step=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Loop delay calculation %%%%%%%%%%%%%


phm = deg2rad(phm_d);

ADC_reso=(ADC_higher_range-ADC_lower_range)/2^ADC_reso_bits;
DPWM_reso=(DPWM_upper_limit-DPWM_lower_limit)/2^DPWM_reso_bits;

D = V_ref/Vin;
Dp = 1-D;

Ts=1/F_SW;
Tsamp=1/N/F_SW;
F_LC = 1/2/pi/sqrt(L*C); F_ESR = 1/2/pi/(ESR_C*C);
a_IRPL = (Vin-V_ref)*D/(2*L*F_SW);
a_OVR = Vin*(1-D)/(16*L*C*(F_SW^2))+2*a_IRPL*ESR_C;

t_const_adc=0;
dts=D*Ts-floor(N*D)/N*Ts;
dpts= Ts - dts;
td=dts;
%%%%%%%%%%%1280%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% State Space Modelling %%%%%%%%%%%%
z = tf('z',Tsamp);
s = tf('s');

A1 = [-(ESR_C+ESR_L)/L -1/L; 1/C 0 ];
A0 = A1;
b1 = [1/L   ESR_C/L; 0 -1/C];
b0 = [0     ESR_C/L;0  -1/C];
c1 = [1 0;ESR_C 1];
c0 = c1;

A1i = A1^-1;
A0i = A0^-1;

Xdown = ((eye(2)-expm(A1*dts)*expm(A0*dpts))^-1)*...
    (-expm(A1*dts)*A0i*(eye(2)-expm(A0*dpts))*b0+...
    -A1i*(eye(2)-expm(A1*dts))*b1)*[Vin;I_Load];
Fdown = (A1-A0)*Xdown + (b1-b0)*[Vin;I_Load];
Phi = expm(A0*(Tsamp-td))*expm(A1*dts)*expm(A0*(td-dts));
gamma =  expm(A0*(Tsamp-td))*Fdown*Tsamp;
delta = c0;

sys = ss(Phi,gamma,delta(2,:),0,Tsamp,'Inputdelay',0);
T = tf(sys);
uncomp = Vin*(1 + s * ESR_C * C)/(1 + s *(ESR_C + ESR_L)*C + s^2 * L * C)*exp(-s*(td));
comp_gain = ADC_reso/DPWM_reso;

%%%%%% L & E optimisation %%%%%%%%%%%%%%%%%%%%%%
if auto == 1
    [ivropts] = LEoptimise(F_SW,L,C,ESR_L,ESR_C,N);
    [~,index] = sortrows([ivropts.eff].');
    ivropts = ivropts(index(end:-1:1));
    clear index
end
%%%%%% PID Compensator realization %%%%%%%%%%%%
flag = 0;
index3 = 1;
if auto == 1
     for (index3 = 1:numel(ivropts))
        flag = 0;
        while(flag == 0)
            F_SW = ivropts(index3).fsw;
            Ts=1/ivropts(index3).fsw;
            Tsamp=1/N/ivropts(index3).fsw;
            L = ivropts(index3).L;
            ESR_L = ivropts(index3).ESR_L;
            F_LC = 1/2/pi/sqrt(L*C); F_ESR = 1/2/pi/(ESR_C*C);
            a_IRPL = (Vin-V_ref)*D/(2*L*F_SW);
            a_OVR = Vin*(1-D)/(16*L*C*(F_SW^2))+2*a_IRPL*ESR_C;
            
            t_const_adc=0;
            dts=D*Ts-floor(N*D)/N*Ts;
            dpts= Ts - dts;
            td=dts;
            %%%%%%%%%%%1280%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%% State Space Modelling %%%%%%%%%%%%
            z = tf('z',Tsamp);
            s = tf('s');
           
            A1 = [-(ESR_C+ESR_L)/L -1/L; 1/C 0 ];
            A0 = A1;
            b1 = [1/L   ESR_C/L; 0 -1/C];
            b0 = [0     ESR_C/L;0  -1/C];
            c1 = [1 0;ESR_C 1];
            c0 = c1;
            
            A1i = A1^-1;
            A0i = A0^-1;
            
            Xdown = ((eye(2)-expm(A1*dts)*expm(A0*dpts))^-1)*...
                (-expm(A1*dts)*A0i*(eye(2)-expm(A0*dpts))*b0+...
                -A1i*(eye(2)-expm(A1*dts))*b1)*[Vin;I_Load];
            Fdown = (A1-A0)*Xdown + (b1-b0)*[Vin;I_Load];
            Phi = expm(A0*(Tsamp-td))*expm(A1*dts)*expm(A0*(td-dts));
            gamma =  expm(A0*(Tsamp-td))*Fdown*Tsamp;
            delta = c0;
            
            sys = ss(Phi,gamma,delta(2,:),0,Tsamp,'Inputdelay',0);
            T = tf(sys);
%             uncomp = Vin*(1 + s * ESR_C * C)/(1 + s *(ESR_C + ESR_L)*C + s^2 * L * C)*exp(-s*(td));
            maxphimd = -1;
            Fc = 0.6*F_SW;
            while(maxphimd<60 & Fc > 0.05*F_SW)
                Wc = 2*pi*Fc;
                [m,p] = bode(uncomp,Wc);
                phimu = 180 + p;
                Wp = 2/Tsamp;
                maxphimd = phimu + 90 - atand(Wc/Wp);
                phm = deg2rad(maxphimd-10);
                Wpd = Wc/(tan(phm - deg2rad(phimu)  + atan(Wc/Wp)));
                Wpi = Wc/50;
                Gpd0 = (1/m)*(sqrt(1 + (Wc/Wp)^2)/(sqrt(1 + (Wc/Wpd)^2)));
                
                Kp = Gpd0*(1+Wpi/Wpd - 2*Wpi/Wp);
                Ki = 2*Gpd0*Wpi/Wp;
                Kd = Gpd0/2*(1-Wpi/Wp)*(Wp/Wpd - 1);
                Gcz = Kp + Ki/(1-z^-1) + Kd*(1-z^-1);
                [Gm,Pm,Wgm,Wpm] = margin(T*Gcz);
                
                
                Wc = 2*pi*Fc;
                [m,p] = bode(T,Wc);
                phimu =  180 + p;
                Wp = 2/Tsamp;
                Wpd = Wc/(tan(phm - deg2rad(phimu)  + atan(Wc/Wp)));
                Wpi = Wc/20;
                Gpd0 = (1/m)*(sqrt(1 + (Wc/Wp)^2)/(sqrt(1 + (Wc/Wpd)^2)));
                maxphimd = phimu + 90 - atand(Wc/Wp);
                Kp = Gpd0*(1+Wpi/Wpd - 2*Wpi/Wp);
                Ki = 2*Gpd0*Wpi/Wp;
                Kd = Gpd0/2*(1-Wpi/Wp)*(Wp/Wpd - 1);
                Gcz = Kp + Ki/(1-z^-1) + Kd*(1-z^-1);
                
                
                b1_val = (Gpd0/2) * ( 1 + Wpi/Wpd + Wp/Wpd + Wpi/Wp);
                b2_val = Gpd0 * (Wpi/Wp - Wp/Wpd);
                b3_val = (Gpd0/2) * (1 - Wpi/Wp)*(Wp/Wpd - 1);
                if Pm > 30
                     simout(Vin,ADC_higher_range,ADC_lower_range,ADC_reso,comp_gain,Tsamp,wk0_bits,fb_wk1_bits,fb_wk2_bits,fb_wk3_bits,...
        b1_val,b2_val,b3_val,Kp,Ki,Kd,DPWM_upper_limit,DPWM_lower_limit,DPWM_reso,F_SW,L,C,ESR_L,ESR_C,load_step_time,load_step,...
        I_load_init,ref_step,vref,sim_time,sim_step_min,sim_step_max);
                    [ivropts(index3).Tset ivropts(index3).Vdroop ivropts(index3).vrt] = calctime(load_step_time,vref);
                    flag = 1;
                end
                Fc = 0.9*Fc;
            end
        end
    end
else
    Wc = 2*pi*Fc;
    [m,p] = bode(T,Wc);
    phimu = 180 + p;
    Wp = 2/Tsamp;
    Wpd = Wc/(tan(phm - deg2rad(phimu)  + atan(Wc/Wp)));
    Wpi = Wc/20;
    Gpd0 = (1/m)*(sqrt(1 + (Wc/Wp)^2)/(sqrt(1 + (Wc/Wpd)^2)));
    maxphimd = phimu + 90 - atand(Wc/Wp);
    Kp = Gpd0*(1+Wpi/Wpd - 2*Wpi/Wp);
    Ki = 2*Gpd0*Wpi/Wp;
    Kd = Gpd0/2*(1-Wpi/Wp)*(Wp/Wpd - 1);
    Gcz = Kp + Ki/(1-z^-1) + Kd*(1-z^-1);
    
    
    b1_val = (Gpd0/2) * ( 1 + Wpi/Wpd + Wp/Wpd + Wpi/Wp);
    b2_val = Gpd0 * (Wpi/Wp - Wp/Wpd);
    b3_val = (Gpd0/2) * (1 - Wpi/Wp)*(Wp/Wpd - 1);
    
end

 Gczs = d2c(Gcz,'tustin');
 loop = uncomp*Gczs;
 Zol = ESR_L * (1 + s * ESR_C * C) * (1 + s*L/ESR_L)/(1 + s *(ESR_C + ESR_L)*C + s^2 * L * C);
 Zocl = Zol/(1+loop);


% b1_val = 	
% fb_wk2_bits		9
% fb_wk3_bits		8
% opts = bodeoptions;
% opts.Magunits = 'abs';
% opts.Magscale = 'log';
% opts.PhaseVisible = 'off';
% opts.FreqUnits = 'Hz';
% 
% opst.Xlim = [2*pi*1e5 2*pi*1e11]
% bode(Zol,opts);
% hold on
% bode(Zocl,opts);
% %
%  freq = logsivropts(index3).Vdrooppace(-3,1,200)*1e9; %THIS IS THE FREQUENCY RANGE THE PDN IMPEDANCE IS CALCULATED
%  [mag phase] = bode(Zocl,2*pi*freq);
%  mag = squeeze(mag);
%  phase = squeeze(phase);
%  Z_VRM = mag.*exp(j.*phase);
% plot(2*pi*freq,abs(Z_VRM));
if (auto ==1)
  save('ivropts.mat','ivropts')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Programming Simulink Model %%%%%%%%%%%
if (simu == 1)
    simout(Vin,ADC_higher_range,ADC_lower_range,ADC_reso,comp_gain,Tsamp,wk0_bits,fb_wk1_bits,fb_wk2_bits,fb_wk3_bits,...
        b1_val,b2_val,b3_val,Kp,Ki,Kd,DPWM_upper_limit,DPWM_lower_limit,DPWM_reso,F_SW,L,C,ESR_L,ESR_C,load_step_time,load_step,...
        I_load_init,ref_step,vref,sim_time,sim_step_min,sim_step_max);
end

