%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Arvind Singh, PhD student in GREEN laboratory 
%Georgia Institue of Technology.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Eff, I_Load, P_TOTAL, V_Ripple]=computeEfficiency(I_total, FSW,L,C,ESR_L,ESR_C, numPhases, rTableIndex,includeFreqTable, VIN, tech, ESR_PCB_PCKG, ESR_Interp, G_Interp, off_chip);
R_m=csvread('table_w_perm_variation.csv');
r=R_m(:,rTableIndex);
f=R_m(:,1);
load_scale=2; %1 - full load opt, 2 - half-load opt
sw_loss_scale=2;
res_loss_scale=2;
if(nargin<8)
    disp('Incorrect Number of Arguments');
    disp('Correct Usage: [Efficiency, I_Load]=computeEfficiency(FSW,L,C,ESR_L,ESR_C, numPhases, rTableIndex, includeFreqTable)');
    disp('includeFreqTable = Y/N');
    disp('rTableIndex=1 -> 6');
else
%     I_total=1.5;
    L_phi=L/numPhases; C=C; ESR_C=ESR_C/numPhases; ESR_L=ESR_L/numPhases; FSW=FSW;
%     Istart=I_total/numPhases; Imax=I_total/numPhases; Istep=I_total;
    Istart=50e-3/numPhases; Imax=I_total/numPhases; Istep=25e-3;
    Iopt = I_total/numPhases/load_scale; % 50% of peak load
    sizeArray=round(((Imax-Istart)/Istep)+1);
    Efficiency=zeros(sizeArray,1);
    V_REF = 1; VRAMP=2;
    RL_phi = 0.01;N = 1;L = L_phi;
    L_BW = 3e-11;ESR_L_BW = ESR_L;C_DECAP = 1e-10;ESR_C_DECAP = ESR_C;comp_gain=50;
    F_Reso = 1/2/pi/sqrt(L*C);
    FO = F_Reso*5;disp(FO/1e6);theta=70;ctr_delay_P = 0e-9;
    tSW=1/FSW;
    
    % stacked topology with 1.8V devices
    if tech=='28nm'
        tech_r_scale=1;
        tech_c_scale=1;
        Cgsn_f=1.011e-12/tech_c_scale; %per mm
        Cgdn_f=0.2443e-12/tech_c_scale; %per mm
        Cgsp_f=0.7968e-12/tech_c_scale; %per mm
        Cgdp_f=0.6735e-12/tech_c_scale; %per mm
        
        r_scale=1;
        Rp_f = 2.147/r_scale/tech_r_scale; % for 1mm
        Rn_f = 2.65/r_scale/tech_r_scale; % for 1mm
    elseif tech=='65nm'
        tech_r_scale=1;
        tech_c_scale=1;
        Cgsn_f=1.1178e-12/tech_c_scale; %per mm
        Cgdn_f=0.3276e-12/tech_c_scale; %per mm
        Cgsp_f=0.9508e-12/tech_c_scale; %per mm
        Cgdp_f=1.0444e-12/tech_c_scale; %per mm
        
        r_scale=1;
        Rp_f = 3.4395/r_scale/tech_r_scale; % for 1mm
        Rn_f = 4.445/r_scale/tech_r_scale; % for 1mm
    end
       
    D_CCM = V_REF/VIN;
    I_Ripple=(VIN-V_REF)*D_CCM/2/L/FSW; % amplitude not p2p
    % disp(['Iopt=' num2str(Iopt)]);
    
    %[w_P, w_N] = OptMOS(VIN, V_REF, Iopt/2, L, 1/FSW);
%     PG_MP_factor = ((Cgsp_f + Cgdp_f)/2)*3*(VIN/2)^2/tSW;
%     PG_MN_factor = ((Cgsn_f + Cgdn_f)/2)*3*(VIN/2)^2/tSW;
    PG_MP_factor = ((2*Cgsp_f + 3*Cgdp_f))*(VIN/2)^2/tSW/sw_loss_scale;
    PG_MN_factor = ((2*Cgsn_f + 3*Cgdn_f))*(VIN/2)^2/tSW/sw_loss_scale;
    PR_MP_factor =(Iopt^2+(I_Ripple^2/3))*D_CCM*Rp_f;
    PR_MN_factor =(Iopt^2+(I_Ripple^2/3))*(1-D_CCM)*Rn_f;
    
    w_P=sqrt(PR_MP_factor/PG_MP_factor);
    w_N=sqrt(PR_MN_factor/PG_MN_factor);
    
    % ratio=2;
    % w_P=w_N*ratio;
    scale=1;
    w_P=scale*w_P;
    w_N=scale*w_N;
%     w_P=169.5971; w_N=347.9839; % off-chip; optimized for VOUT=1V and then fixed for all votlages
%     w_P=21.5057; w_N=44.1261; % on-chip; optimized for VOUT=1V and then fixed for all votlages
%     w_P=62.5; w_N=128.3; % on-chip; optimized for VOUT=0.8 and then fixed
%     for all voltages. for 5X higher L, C (L=120n, C=110n). FSW=15M.
%     w_P=48.1; w_N=98.7; % on-chip; optimized for VOUT=0.8 and then fixed for all voltages. 
    % for 5X higher L, C (L=120n, C=110n) FSW=25M.

    Rds_on_P = Rp_f/w_P;
    Rds_on_N = Rn_f/w_N;
    Cpfet = (Cgsp_f+Cgdp_f)*w_P;
    Cnfet = (Cgsn_f+Cgdn_f)*w_N;
    
%     disp(['w_P is ' num2str(w_P) ' w_N is ' num2str(w_N)]);
%     disp(['Rds_on_P is ' num2str(Rds_on_P) ' Rds_on_N is ' num2str(Rds_on_N)]);
%     disp(['Cnfet is ' num2str(Cnfet) ' Cpfet is ' num2str(Cpfet)]);
%     disp(['D_CCM is ' num2str(D_CCM)]);
%     disp(['I_Ripple is ' num2str(I_Ripple)]);
    for I_Load=Istart:Istep:Imax;
        index=round(((I_Load-Istart)/Istep)+1);
%         disp(['I_Load=' num2str(I_Load)]);
        if(includeFreqTable=='Y')
            % the below is from Sebastian's simulations. a factor of 2 comes from +/- freq. components
            P_INDUCTOR_AC=2*I_Ripple^2*(0.405^2*(interp1(f,r,FSW)) + 0.045^2*(interp1(f,r,3*FSW)) + 0.016^2*(interp1(f,r,5*FSW)));
        else
            P_INDUCTOR_AC=2*I_Ripple^2*(0.405^2*(interp1(f,r,0)) + 0.045^2*(interp1(f,r,0)) + 0.016^2*(interp1(f,r,0)));
        end
        
        % P_INDUCTOR_DC=I_Load^2*interp1(f,r,0);
        P_INDUCTOR_DC=I_Load^2*ESR_L;
        P_FET_RESISTIVE=res_loss_scale*((I_Load^2+(I_Ripple^2)/3)*D_CCM*Rp_f/w_P + (I_Load^2+(I_Ripple^2)/3)*(1-D_CCM)*Rn_f/w_N); % factor of 2 for stacked topology
        P_FET_SWITCHING=sw_loss_scale*(PG_MP_factor*w_P+ PG_MN_factor*w_N);
        P_CAP_LOSS_AC=I_Ripple^2/3*ESR_C;    
        if off_chip=='Y'
        P_PCB_PCKG=I_Load^2*ESR_PCB_PCKG;%+V_REF^2*G_PCB_PCKG;
    else
        P_PCB_PCKG=I_Load^2*ESR_PCB_PCKG*(V_REF/VIN)^2;%+VIN^2*G_PCB_PCKG;
    end
        P_Interposer=I_Load^2*ESR_Interp;%+V_REF^2*G_Interp;

%         disp(['FETs RESISTIVE is ' num2str(P_FET_RESISTIVE)]);
%         disp(['FETs SWITCHING is ' num2str(P_FET_SWITCHING)]);
%         disp(['INDUCTOR DC is ' num2str(P_INDUCTOR_DC)]);
%         disp(['INDUCTOR AC is ' num2str(P_INDUCTOR_AC)]);
%         disp(['PCB+PCKG is ' num2str(P_PCB_PCKG)]);
%         disp(['Interposer is ' num2str(P_Interposer)]);
%         disp(['ESR CAP is ' num2str(P_CAP_LOSS_AC)]);
        
        P_TOTAL = P_FET_RESISTIVE + P_FET_SWITCHING + P_INDUCTOR_DC + P_INDUCTOR_AC + P_PCB_PCKG + P_Interposer + P_CAP_LOSS_AC;
        
        Efficiency(index,1)=(V_REF*I_Load)/(P_TOTAL+V_REF*I_Load);
    end
    I_Load=Istart:Istep:Imax;
    [y1, index]=max(Efficiency*100);
%     I_opt=((index-1)*Istep)+Istart;
    I_opt=I_total/numPhases/load_scale;
    I_Load=I_opt;
%     I_Load=I_total/numPhases/2;
    
    if(includeFreqTable=='Y')
        % the below is from Sebastian's simulations. a factor of 2 comes from +/- freq. components
        P_INDUCTOR_AC=2*I_Ripple^2*(0.405^2*(interp1(f,r,FSW)) + 0.045^2*(interp1(f,r,3*FSW)) + 0.016^2*(interp1(f,r,5*FSW)));
    else
        P_INDUCTOR_AC=2*I_Ripple^2*(0.405^2*(interp1(f,r,0)) + 0.045^2*(interp1(f,r,0)) + 0.016^2*(interp1(f,r,0)));
    end
    
    % P_INDUCTOR_DC=I_Load^2*interp1(f,r,0);
    P_INDUCTOR_DC=I_Load^2*ESR_L;
    P_Interposer=I_Load^2*ESR_Interp; %+V_REF^2*G_Interp;
    P_FET_RESISTIVE=res_loss_scale*(((((I_Load^2)+((I_Ripple^2)/3))*D_CCM*Rp_f)/w_P) + (((((I_Load^2)+(I_Ripple^2)/3))*(1-D_CCM)*Rn_f)/w_N)); % factor of 2 for stacked topology
    P_FET_SWITCHING=sw_loss_scale*(PG_MP_factor*w_P+ PG_MN_factor*w_N);
    P_CAP_LOSS_AC=I_Ripple^2/3*ESR_C;
    if off_chip=='Y'
        P_PCB_PCKG=I_Load^2*ESR_PCB_PCKG;%+V_REF^2*G_PCB_PCKG;
    else
        P_PCB_PCKG=I_Load^2*ESR_PCB_PCKG*(V_REF/VIN)^2;%+VIN^2*G_PCB_PCKG;
    end
%     V_Ripple=(1/2)*(I_Ripple)*(1/(2*FSW))/C/numPhases; % doesn't model
%     impact of ESR_C
    V_Ripple=(1/2)*(I_Ripple)*(1/(2*FSW))/C/numPhases + 2*I_Ripple*ESR_C/numPhases; % linear model - p2p
%     V_Ripple=sqrt(((1/2)*(I_Ripple)*(1/(2*FSW))/C/numPhases)^2 +   (2*I_Ripple*ESR_C)^2/numPhases); % rms model: p2p
    
    disp(['%%%%% Optimal Point is ' num2str(I_opt) 'A %%%%%']);
    disp(['w_P is ' num2str(w_P) ' w_N is ' num2str(w_N)]);
    disp(['Rds_on_P is ' num2str(Rds_on_P) ' Rds_on_N is ' num2str(Rds_on_N)]);
    disp(['Cnfet is ' num2str(Cnfet) ' Cpfet is ' num2str(Cpfet)]);
    disp(['D_CCM is ' num2str(D_CCM)]);
    disp(['I_Ripple pk-pk is ' num2str(2*I_Ripple)]);
    disp(['V_Ripple pk-pk is ' num2str(V_Ripple)]);
    disp(['FETs RESISTIVE is ' num2str(P_FET_RESISTIVE)]);
    disp(['FETs SWITCHING is ' num2str(P_FET_SWITCHING)]);
    disp(['INDUCTOR DC is ' num2str(P_INDUCTOR_DC)]);
    disp(['INDUCTOR AC is ' num2str(P_INDUCTOR_AC)]);
    disp(['PCB+PCKG is ' num2str(P_PCB_PCKG)]);
    disp(['ESR CAP is ' num2str(P_CAP_LOSS_AC)]);
    disp(['Interposer is ' num2str(P_Interposer)]);
        
    P_TOTAL = P_FET_RESISTIVE + P_FET_SWITCHING + P_INDUCTOR_DC + P_INDUCTOR_AC + P_PCB_PCKG + P_Interposer + P_CAP_LOSS_AC;
    
    Efficiency(index,1)=(V_REF*I_Load)/(P_TOTAL+V_REF*I_Load);
    disp(['Power Loss ' num2str(P_TOTAL)]);
    disp(['Efficiency ' num2str((V_REF*I_Load)/(P_TOTAL+V_REF*I_Load))]);
    Eff = (V_REF*I_Load)/(P_TOTAL+V_REF*I_Load);
    I_Load=Istart:Istep:Imax;
%     figure(20);plot(I_Load,Efficiency,'Linewidth',2);
end
end
