%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Nihar Dasari, PhD student in GREEN laboratory 
%Georgia Institue of Technology.
%Function to simulate IVR in Simulink
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time] =  simout(Vin,ADC_higher_range,ADC_lower_range,ADC_reso,comp_gain,Tsamp,wk0_bits,fb_wk1_bits,fb_wk2_bits,fb_wk3_bits,...
    b1_val,b2_val,b3_val,Kp,Ki,Kd,DPWM_upper_limit,DPWM_lower_limit,DPWM_reso,F_SW,L,C,ESR_L,ESR_C,load_step_time,load_step,...
    I_load_init,ref_step,vref,sim_time,sim_step_min,sim_step_max);
   

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%%%%%%% Programming Simulink Model %%%%%%%%%%%
    file_str='IVR';
    open_system(file_str);
    
    % Err Out ADC
    set_param(strcat(file_str,'/Buck/ZOH'),'SampleTime',num2str(Tsamp));
    % set_param(strcat(file_str,'/Buck/ADC_out_delay'),'SampleTime',num2str(Tsamp));
    set_param(strcat(file_str,'/Buck/ADC_out'),'QuantizationInterval',num2str(ADC_reso));
    set_param(strcat(file_str,'/Buck/lim_ADC_out'),'UpperLimit',num2str(ADC_higher_range));
    set_param(strcat(file_str,'/Buck/lim_ADC_out'),'LowerLimit',num2str(ADC_lower_range));
    set_param(strcat(file_str,'/Buck/ADC_to_bin'),'Gain',num2str(1/ADC_reso));
    
    % Compensator
    set_param(strcat(file_str,'/Buck/Comp_Gain'),'Gain',num2str(comp_gain));
    
    % PID direct form %
%     set_param(strcat(file_str,'/Buck/PID_direct/Delay_first'),'SampleTime',num2str(Tsamp));
%     set_param(strcat(file_str,'/Buck/PID_direct/Delay_second'),'SampleTime',num2str(Tsamp));
%     
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_wk0'),'UpperLimit',num2str(2^wk0_bits-1));
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_wk0'),'LowerLimit',num2str(-2^wk0_bits));
%     
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_ffwk1'),'UpperLimit',num2str(2^fb_wk1_bits-1));
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_ffwk1'),'LowerLimit',num2str(-2^fb_wk1_bits));
%     
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_ffwk2'),'UpperLimit',num2str(2^fb_wk2_bits-1));
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_ffwk2'),'LowerLimit',num2str(-2^fb_wk2_bits));
%     
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_ffwk3'),'UpperLimit',num2str(2^fb_wk3_bits-1));
%     set_param(strcat(file_str,'/Buck/PID_direct/lim_ffwk3'),'LowerLimit',num2str(-2^fb_wk3_bits));
%     
%     set_param(strcat(file_str,'/Buck/PID_direct/b1'),'Gain',num2str(b1_val));
%     set_param(strcat(file_str,'/Buck/PID_direct/b2'),'Gain',num2str(b2_val));
%     set_param(strcat(file_str,'/Buck/PID_direct/b3'),'Gain',num2str(b3_val));
%     
    % Backward Euler PID parallel form %
    Kp_1=Kp;
    Ki_1=Ki/Tsamp;
    Kd_1=Kd*Tsamp;
    
    set_param(strcat(file_str,'/Buck/PID_parallel'),'SampleTime',num2str(Tsamp));
    set_param(strcat(file_str,'/Buck/PID_parallel'),'P',num2str(Kp_1));
    set_param(strcat(file_str,'/Buck/PID_parallel'),'I',num2str(Ki_1));
    set_param(strcat(file_str,'/Buck/PID_parallel'),'D',num2str(Kd_1));
    
    % DPWM
    % set_param(strcat(file_str,'/Buck/DPWM'),'SampleTime',num2str(Ts));
    % Total DPWM range is 2 (-1 to 1) and DPWM input bits should cover the
    % entire range
    set_param(strcat(file_str,'/Buck/lim_dpwm_out'),'UpperLimit',num2str(DPWM_upper_limit));
    set_param(strcat(file_str,'/Buck/lim_dpwm_out'),'LowerLimit',num2str(DPWM_lower_limit));
    % set_param(strcat(file_str,'/Buck-1/bin_quant'),'SampleTime',num2str(Tsamp));
    set_param(strcat(file_str,'/Buck/bin_to_DPWM'),'Gain',num2str(DPWM_reso));
    set_param(strcat(file_str,'/Buck/PWM/Sawtooth'),'Freq',num2str(F_SW));
    
    
    % Power stage  < numel(ivropts)
    set_param(strcat(file_str,'/Input_voltage'),'Value',num2str(Vin));
    set_param(strcat(file_str,'/Buck/buck_converter/1_by_L_filt'),'Gain',num2str(1/L));
    set_param(strcat(file_str,'/Buck/buck_converter/1_by_C_filt'),'Gain',num2str(1/C));
    set_param(strcat(file_str,'/Buck/buck_converter/RL_filt'),'Gain',num2str(ESR_L));
    set_param(strcat(file_str,'/Buck/buck_converter/Resr_filt'),'Gain',num2str(ESR_C));
%     set_param(strcat(file_str,'/Buck/buck_converter/1_by_L_BW'),'Gain',num2str(1/L_BW));
%     set_param(strcat(file_str,'/Buck/buck_converter/1_by_C_PAD'),'Gain',num2str(1/C_DECAP));
%     set_param(strcat(file_str,'/Buck/buck_converter/Resr_PAD'),'Gain',num2str(ESR_C_DECAP));
%     set_param(strcat(file_str,'/Buck/buck_converter/RL_BW'),'Gain',num2str(ESR_L_BW));
    
    % Trasient events
    set_param(strcat(file_str,'/Output_Load/Delay_Load'),'DelayTime',num2str(load_step_time));
    set_param(strcat(file_str,'/Output_Load/Pulse Generator'),'Amplitude',num2str(load_step));
    set_param(strcat(file_str,'/Output_Load/I_ref'),'Value',num2str(I_load_init));
    
    set_param(strcat(file_str,'/Buck/Pulse Generator2'),'Amplitude',num2str(ref_step));
    set_param(strcat(file_str,'/Vref'),'Value',num2str(vref));
    
    
    set_param(file_str,'StopTime',num2str(sim_time),'MinStep',num2str(sim_step_min),'MaxStep',...
        num2str(sim_step_max),'SimulationMode','normal');
   sim(file_str,[]);
    save_system(file_str);
%     close_system(file_str)  
end