%% Generate virtual cell population
clear all; clc; close all;

%parameters for sampling
tempconvt_nominal = [0.13, 0.65, 0.65,... %kin0, kin1, kout1
                     8, 20, 0.45]; %kin2, Kin2, kout2 
hsentrain_nominal = [27.1731, 3.6760, 39.9689, 1.0723,... %ks,Ks,khsf1,Khsf1
                     20.8923, 0.3343, 20.6516, 8.7180, 2.3686]; %HSF1tot, vact, vina, vbhsp, vdhsp
pcg_nominal = [9, 1, 0.56, 0.12, ... %v1b, k1b, k1i, k1d
               0.3, 0.05, 0.24, ... %k2b, k2d, k2t 
               0.02, 0.12, ... %k3t, k3d
               3.6, 2.16, 0.75, ... %v4b, k4b, k4d
               0.24, 0.06, 0.45, ... %k5b, k5d, k5t
               0.06, 0.12, 0.09, ... %k6t, k6d, k6a
               0.003, 0.09]; %k7a, k7d
var_center = [tempconvt_nominal,hsentrain_nominal,pcg_nominal];
var_min = var_center.*(1-0.03); %range for the sampling +/-3%
var_max = var_center.*(1+0.03);

%construct d-D Sobol quasi-random point set for sampling
ps = sobolset(length(var_center)); 
npop = 1000; %cell population size
pop = [];
for i = 1:npop
    pop0 = var_min+(var_max-var_min).*ps(i,:);
    pop = [pop;pop0];
end

writematrix(pop,'pop_ALL_samp3%.xlsx','WriteMode','overwrite');


%% Fig.2: temporal dynamics, phase, phase relationship under W12/C12 37±1.5°C temperature rhythm (single cell, nominal parameter set)
clear all; close all; clc; 

tRange = [0 4800+24*200];
Y0 = ones(1,11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%(human core body) temperature fluctuation (same to Fig.1 (Saini et al, 2012))
pt = 24;
pw = 12;
Tempint = 3; %°C
Tcoreavg = 37; %°C

temppar = [27.1731  3.6760  39.9689  1.0723  20.8923  0.3343  20.6516  8.7180  2.3686]; 
entrainparam = temppar(1:4);
hs = temppar(5:end);
pcg = [9, 1, 0.56, 0.12, ... %v1b, k1b, k1i, k1d
       0.3, 0.05, 0.24, ... %k2b, k2d, k2t 
       0.02, 0.12, ... %k3t, k3d
       3.6, 2.16, 0.75, ... %v4b, k4b, k4d
       0.24, 0.06, 0.45, ... %k5b, k5d, k5t
       0.06, 0.12, 0.09, ... %k6t, k6d, k6a
       0.003, 0.09]; %k7a, k7d

[t,Y] = ode45(@(t,Y)Temp_HSF1_PCG(t,Y,pw,pt,Tempint,Tcoreavg,pcg,hs,entrainparam),tRange,Y0,options);
mpercry = Y(:,1); percry = Y(:,2); nucpercry = Y(:,3); 
mbmal1 = Y(:,4); bmal1 = Y(:,5); nucbmal1 = Y(:,6);
clockbmal1 = Y(:,7);
acthsf1 = Y(:,8); indhsp = Y(:,9);
ts1 = Y(:,10); ts2 = Y(:,11); 

%temporal dynamics
subplot(2,2,[1,2])
line_color = {[52 57 60]./255, [219 73 76]./255, [235 164 50]./255, [134 170 109]./255, [70 181 211]./255, [102 103 171]./255};
yyaxis left
zz0 = Tcoreavg-Tempint/2;
zz1 = Tcoreavg+Tempint/2;
x1 = [0 24;    12 36;   12 36;   0 24];
z1 = [zz0 zz0; zz0 zz0; zz1 zz1; zz1 zz1];
f1 = fill(x1,z1,[130 19 45]./255,'edgecolor','none'); alpha(0.15); hold on

x2 = [12 36; 24 48; 24 48; 12 36];
z2 = [zz0 zz0; zz0 zz0; zz1 zz1; zz1 zz1];
f2 = fill(x2,z2,[68 89 168]./255,'edgecolor','none'); alpha(0.15); hold on
ylabel('Temperature (°C)'); 
set(gca,'YColor',line_color{2});
        
yyaxis right
pl(1) = plot(t-4800,ts2,'-','Color',line_color{3},'LineWidth',2);
pl(2) = plot(t-4800,acthsf1,'--','Color',line_color{3},'LineWidth',2);
pl(3) = plot(t-4800,mpercry,'-','Color',line_color{4},'LineWidth',2); hold on
pl(4) = plot(t-4800,percry,'--','Color',line_color{4},'LineWidth',2); hold on
pl(5) = plot(t-4800,mbmal1,'-','Color',line_color{5},'LineWidth',2); hold on
pl(6) = plot(t-4800,bmal1,'--','Color',line_color{5},'LineWidth',2); hold on
xlim([0 48]); xticks(0:6:48); xlabel('Time (h)'); ylabel('Concentration (a.u.)'); 
legend([f1(1),f2(1),pl],{'warm period','cold period','TS_2','activated HSF1','Per/Cry mRNA','PER/CRY','Bmal1 mRNA','BMAL1'},'Location','north','Orientation','horizontal')
legend('boxoff')
set(gca,'YColor',line_color{1},'FontSize',16);


%EXTRA: period
idx_pd = find(t>=4800); 
t_pd = t(idx_pd);
pd0 = [];
for i = 1:length(Y0)
    compn = Y(idx_pd,i);
    idmax_pd = find(diff(sign(diff(compn)))<0)+1;
    pd0(end+1) = mean(diff(t_pd(idmax_pd)));
end
compn_name = ["Entrainer period (Calibration)";"Per/Cry mRNA";"PER/CRY";"nucleus PER/CRY";"Bmal1 mRNA";"BMAL1";"nucleus BMAL1";"CLOCK/BMAL1";"actHSF1";"indHSP";"TS1";"TS2"];
pd_val = [pt;pd0'];
unit = ['hr';'hr';'hr';'hr';'hr';'hr';'hr';'hr';'hr';'hr';'hr';'hr'];
pd_tbl = table(compn_name,pd_val,unit)


%phase calibration of PER/CRY(pred) with PER2(exp), phase relationship between mPer/Cry and mBmal1
index = find(t>=4800 & t<4800+49);
t0 = t(index);
mpercry0 = mpercry(index); 
percry0 = percry(index); 
mbmal10 = mbmal1(index);

indmax_mpercry = []; indmax_percry = []; indmax_mbmal1 = [];
for l = 2:length(t0)-1
    if mpercry0(l-1) < mpercry0(l) & mpercry0(l+1) < mpercry0(l) 
        indmax_mpercry(end+1) = l;
    end
    if percry0(l-1) < percry0(l) & percry0(l+1) < percry0(l)
        indmax_percry(end+1) = l;
    end
    if mbmal10(l-1) < mbmal10(l) & mbmal10(l+1) < mbmal10(l)
        indmax_mbmal1(end+1) = l;
    end
end    
pkt_percry = mod(t0(indmax_percry(1)),24);
pkt_cr_mpercry = mod(t0(indmax_mpercry(1))-4812,24); %include a "correction" for aligning with the temp cycle used in Fig.1 (Saini et al, 2012), which is starting at the ~24+12 hour
pkt_cr_mbmal1 = mod(t0(indmax_mbmal1(1))-4812,24);

subplot(2,2,3)
data_exp_cosinor = xlsread('data_exp_cosinor.xlsx'); %post-cosinor data of Fig.1C (Saini et al, 2012)
data_exp_cosinor = data_exp_cosinor';
t_exp = data_exp_cosinor(1,:)+4800; 
PER_exp = data_exp_cosinor(2,:);

pkt_expPER2 = t_exp(find(PER_exp == max(PER_exp)));
theta_expmodel = [mod(pkt_expPER2,24)/24*(2*pi) pkt_percry/24*(2*pi)]; %exp PER2, pred PER/CRY
rho_percry = [1.35,1.45];
polarscatter(theta_expmodel(1),rho_percry(1),100,'ro','filled','MarkerFaceAlpha',.75,'MarkerEdgeColor','k'); hold on 
polarscatter(theta_expmodel(2),rho_percry(2),100,'ko','filled','MarkerFaceAlpha',.75,'MarkerEdgeColor','k');
title('Phase Calibration');
legend('Experimental PER2, mPer2, mBmal1','Predicted PER/CRY, mPer/Cry, mBmal1','Location','southeastoutside');
set(gca,'RAxisLocation',0,'RTick',[0 1.5],'RTickLabel',{'','PER'},'Rcolor','k','GridAlpha',1,'ThetaZeroLocation','top','ThetaTick',[0 45 90 135 180 225 270 315],'ThetaTickLabel',{'0','3','6','9','12','15','18','21'},'ThetaDir','clockwise','ThetaGrid','off','Color',[0.95 0.95 0.95],'LineWidth',1,'FontSize',16);
        
subplot(2,2,4)
pkt_cir_mpercry_exp = 23; pkt_cir_mbmal1_exp = 12; %Fig.1D (Saini et al, 2012)
theta_exp = [pkt_cir_mpercry_exp/24*(2*pi) pkt_cir_mbmal1_exp/24*(2*pi)];
theta_model = [pkt_cr_mpercry/24*(2*pi) pkt_cr_mbmal1/24*(2*pi)];
rho_mper = [0.75,0.85,0.95];
rho_mbmal1 = [1.75,1.85,1.95];
polarscatter(theta_exp,[rho_mper(1),rho_mbmal1(1)],100,'ro','filled','MarkerFaceAlpha',.75,'MarkerEdgeColor','k'); hold on
polarscatter(theta_model,[rho_mper(2),rho_mbmal1(2)],100,'ko','filled','MarkerFaceAlpha',.75,'MarkerEdgeColor','k');
title('Phase Relationship Calibration');
set(gca,'RAxisLocation',0,'RTick',[0 1 2],'RTickLabel',{'','mPer','mBmal1'},'Rcolor','k','GridAlpha',1,'ThetaZeroLocation','top','ThetaTick',[0 45 90 135 180 225 270 315],'ThetaTickLabel',{'0','3','6','9','12','15','18','21'},'ThetaDir','clockwise','ThetaGrid','off','Color',[0.95 0.95 0.95],'LineWidth',1,'FontSize',16);


%% Fig.3: entrainment responses of cell clock population to normal/nominal and inversed temperature rhythms
clear all; close all; clc;

pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop = pop0;

tRange = 0:1e-1:1200+24*135; %[0 1200+24*135];
Y0 = ones(1*size(pop,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%characteristics of temperature rhythms
pt = 24;
pw = 12;
Tcoreavg = 37;
Tempint = 3;
p_test = [1 -1]; %temperature schedules including either normal/nominal rhythm or inversed rhythm

phasedistri_mpercry = []; perioddistri_mpercry = [];
Rsyn = []; %1st column - Per/Cry mRNA; 2nd column - Bmal1 mRNA

line_typ = {'k-','r:'}; line_width = [1,0.5];
marker_typ = {'ko--','ro--'};

pn = 1
for i = 1:length(p_test)
    
    p = p_test(i); 
    tic;
    [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync(t,Y,pw,pt,Tempint,Tcoreavg,pop,p),tRange,Y0,options);
    toc
    
    if i == 1
        save('ode45Output-normal.mat','t','Y');
    else
        save('ode45Output-inversed.mat','t','Y');
    end
    
    YY = reshape(Y,[],11,size(pop,1));
    mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
    mpercry_ens = mean(mpercry_pop,2);
   
    figure(1) %Fig.3A
    subplot(2,3,[1,2,3]); plot(t-1200,mpercry_ens,line_typ{i},'LineWidth',line_width(i)); hold on; 
    
    %phase and period (distribution), R synchronization index
    tpoint = 0:3:(t(end)-1200)/24-3;
    t_lb_set = 1200+24*tpoint;
    t_ub_set = 1200+24*(tpoint+3);
    pn1 = 1
    for tn = 1:length(t_lb_set)
        t_lb = t_lb_set(tn); t_ub= t_ub_set(tn);
        index = find(t>=t_lb & t<t_ub);
        t0 = t(index);
        mpercry_pop0 = mpercry_pop(index,:);
        
        phase_mpercry = []; period_mpercry = [];
        pn2 = 1
        for n = 1:size(pop,1)  
            indmax_mpercry = find(diff(sign(diff(mpercry_pop0(:,n))))<0)+1;
            period0_mpercry = mean(diff(t0(indmax_mpercry)));
            period_mpercry = [period_mpercry,period0_mpercry];
            
            indphase_mpercry = indmax_mpercry(1);
            tphase_mpercry = t0(indphase_mpercry);
            phase0_mpercry = (tphase_mpercry-t_lb)/period0_mpercry*24;
            phase_mpercry = [phase_mpercry,phase0_mpercry];
            
            pn2 = pn2 + 1
        end
        
        if t_lb>=1200+24*10 & t_lb<1200+24*80
            phase_mpercry(phase_mpercry>=21) = 24-phase_mpercry(phase_mpercry>=21);
        end
        
        mn1 = mean(phase_mpercry); 
        stdv1 = std(phase_mpercry); 
        phasedistri_mpercry = [phasedistri_mpercry; mn1,stdv1];
        mn2 = mean(period_mpercry);
        stdv2 = std(period_mpercry);
        perioddistri_mpercry = [perioddistri_mpercry; mn2,stdv2];
        
        
        Rpop_mpercry = mpercry_pop(index,:);
        Rt = t(index)-t_lb;

        y_square_mpercry = Rpop_mpercry.^2;
        y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
        y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
        y_tavg_square_mpercry = y_tavg_mpercry.^2;

        ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
        ypopavg_square_mpercry = ypopavg_mpercry.^2;
        ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
        ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
        ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

        Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
        Rsyn = [Rsyn;Rsyn_mpercry];
        
        pn1 = pn1 + 1
    end
    
    figure(1)
    subplot(2,3,4);
    errorbar((t_lb_set-1200)./24,phasedistri_mpercry(1+(i-1)*length(t_lb_set):end,1),phasedistri_mpercry(1+(i-1)*length(t_lb_set):end,2),marker_typ{i},'LineWidth',1); hold on
    xlim([0 120]); xlabel('Time (days)'); ylabel({('Per/Cry mRNA phase distribution'),('(mean \pm SD)')}); 
    set(gca,'XTick',0:20:120,'XTickLabel',{'0','20','40','60','80','100','120'},'FontSize',14);
    box on; axis square

    subplot(2,3,5);
    errorbar((t_lb_set-1200)./24,perioddistri_mpercry(1+(i-1)*length(t_lb_set):end,1),perioddistri_mpercry(1+(i-1)*length(t_lb_set):end,2),marker_typ{i},'LineWidth',1); hold on
    xlim([0 120]); xlabel('Time (days)'); ylabel({('Per/Cry mRNA period distribution'),('(mean \pm SD)')}); 
    set(gca,'XTick',0:20:120,'XTickLabel',{'0','20','40','60','80','100','120'},'FontSize',14);
    box on; axis square
        
    subplot(2,3,6);
    plot((t_lb_set-1200)./24,Rsyn(1+(i-1)*length(t_lb_set):end),marker_typ{i},'LineWidth',1); hold on
    xlim([0 120]); xlabel('Time (days)'); ylabel('R_{syn} for Per/Cry mRNA'); 
    set(gca,'XTick',0:20:120,'XTickLabel',{'0','20','40','60','80','100','120'},'FontSize',14);
    box on; axis square
        
    pn = pn + 1
end

figure(1) %Fig.3A
subplot(2,3,[1,2,3]);
xlim([0 24*120]); ylim([0 2.2]); xlabel('Time (days)'); ylabel({('Per/Cry mRNA'),('ensemble average level (a.u.)')}); 
set(gca,'XTick',0:24*10:24*120,'XTickLabel',{'0','10','20','30','40','50','60','70','80','90','100','110','120'},'FontSize',16);
legend({'normal T rhythm','reversed T rhythm'},'Orientation','vertical','Location','northeast','FontSize',18);  
legend('boxoff');

%% Fig.4: amplitude-varying temperature rhythms
clear all; close all; clc;

pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop = pop0;

tRange = 0:1e-1:2400+24*50; %[0 2400+24*50];
Y0 = ones(1*size(pop,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%characteristics of temperature rhythms
pt = 24;
pw = 12;
Tcoreavg = 37; %°C
p = 0; %type of temperature schedule

tempamp_test = [0,1:1:5];
tempamp_name = {'k_{hsf1} = 0',['T_{amp} = 0 °C' newline '(constant 37 °C)'],'T_{amp} = 1 °C','T_{amp} = 2 °C','T_{amp} = 3 °C','T_{amp} = 4 °C','T_{amp} = 5 °C'};
line_color = {[52 57 60]./255, [52 57 60]./255, [219 73 76]./255, [235 164 50]./255, [134 170 109]./255, [70 181 211]./255, [102 103 171]./255};
phasedistri_acthsf1 = []; perioddistri_acthsf1 = []; 
phasedistri_mpercry = []; perioddistri_mpercry = [];
Rsyn = []; 

rn = 1
for i = 1:length(tempamp_test)+1

    if i == 1
        pop(:,9) = 0;  %w/o entraining effect
        Tempint = 3; %keep amplitude at nominal value of 3 °C
        line_typ = '--';
    else
        pop = pop0; %w/ entraining effect
        Tempint = tempamp_test(i-1);
        line_typ = '-';
    end

    tic;
    [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync(t,Y,pw,pt,Tempint,Tcoreavg,pop,p),tRange,Y0,options);
    toc
    YY = reshape(Y,[],11,size(pop,1));
    mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
    mpercry_ens = mean(mpercry_pop,2);
    acthsf1_pop = reshape(YY(:,8,:),[],size(pop,1));
    acthsf1_ens = mean(acthsf1_pop,2);
    
    figure(1) %Fig.4A
    plot(t-2400,acthsf1_ens,line_typ,'LineWidth',2,'Color',line_color{i}); hold on; 
    xlim([0 48]); ylim([0.6 1.6]); xlabel('Time (h)'); ylabel({('actHSF1'),('ensemble average level (a.u.)')}); set(gca,'XTick',0:6:48,'FontSize',16);
    box on; axis square
    figure(2) %Fig.4B
    plot(t-2400,mpercry_ens,line_typ,'LineWidth',2,'Color',line_color{i}); hold on; 
    xlim([0 48]); ylim([0 2]); xlabel('Time (h)'); ylabel({('Per/Cry mRNA'),('ensemble average level (a.u.)')}); set(gca,'XTick',0:6:48,'FontSize',16); 
    box on; axis square
    
    t_interv = 2400;
    index = find(t>=t_interv);
    t0 = t(index);
    acthsf1_pop0 = acthsf1_pop(index,:);
    mpercry_pop0 = mpercry_pop(index,:);
    
    phase_mpercry = []; period_mpercry = []; 
    phase_acthsf1 = []; period_acthsf1 = [];
    
    if i == 2
        pn = 1
            for n = 1:size(pop,1)  
                indmax_mpercry = find(diff(sign(diff(mpercry_pop0(:,n))))<0)+1;
                period0_mpercry = mean(diff(t0(indmax_mpercry)));
                period_mpercry = [period_mpercry,period0_mpercry];
                
                indphase_mpercry = indmax_mpercry(1);
                tphase_mpercry = t0(indphase_mpercry);
                phase0_mpercry = mod(tphase_mpercry,24);
                phase_mpercry = [phase_mpercry,phase0_mpercry];
                
                pn = pn + 1
            end
            
            mn1 = mean(phase_mpercry); 
            stdv1 = std(phase_mpercry); 
            phasedistri_mpercry = [phasedistri_mpercry; mn1,stdv1];
            mn2 = mean(period_mpercry);
            stdv2 = std(period_mpercry); 
            perioddistri_mpercry = [perioddistri_mpercry; mn2,stdv2];
            
            
            Rpop_mpercry = mpercry_pop(index,:);
            Rt = t(index)-t_interv;

            y_square_mpercry = Rpop_mpercry.^2;
            y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
            y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
            y_tavg_square_mpercry = y_tavg_mpercry.^2;

            ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
            ypopavg_square_mpercry = ypopavg_mpercry.^2;
            ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
            ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
            ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

            Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
            Rsyn = [Rsyn;nan,Rsyn_mpercry];
    else
        pn = 1
            for n = 1:size(pop,1)  
                indmax_acthsf1 = find(diff(sign(diff(acthsf1_pop0(:,n))))<0)+1;
                indmax_mpercry = find(diff(sign(diff(mpercry_pop0(:,n))))<0)+1;
                period0_acthsf1 = mean(diff(t0(indmax_acthsf1)));
                period_acthsf1 = [period_acthsf1,period0_acthsf1];
                period0_mpercry = mean(diff(t0(indmax_mpercry)));
                period_mpercry = [period_mpercry,period0_mpercry];
                
                indphase_acthsf1 = indmax_acthsf1(1);
                tphase_acthsf1 = t0(indphase_acthsf1);
                phase0_acthsf1 = mod(tphase_acthsf1,24);
                phase_acthsf1 = [phase_acthsf1,phase0_acthsf1];
                indphase_mpercry = indmax_mpercry(1);
                tphase_mpercry = t0(indphase_mpercry);
                phase0_mpercry = mod(tphase_mpercry,24);
                phase_mpercry = [phase_mpercry,phase0_mpercry];
                
                pn = pn + 1
            end
            mn1 = mean(phase_mpercry); %calculate the mean
            stdv1 = std(phase_mpercry); %calculate the standard deviation
            phasedistri_mpercry = [phasedistri_mpercry; mn1,stdv1];
            mn2 = mean(period_mpercry);
            stdv2 = std(period_mpercry); 
            perioddistri_mpercry = [perioddistri_mpercry; mn2,stdv2];
            mn3 = mean(phase_acthsf1); 
            stdv3 = std(phase_acthsf1); 
            phasedistri_acthsf1 = [phasedistri_acthsf1; mn3,stdv3];
            mn4 = mean(period_acthsf1);
            stdv4 = std(period_acthsf1); 
            perioddistri_acthsf1 = [perioddistri_acthsf1; mn4,stdv4];
            
            
            Rpop_acthsf1 = acthsf1_pop(index,:);
            Rpop_mpercry = mpercry_pop(index,:);
            Rt = t(index)-t_interv;

            y_square_acthsf1 = Rpop_acthsf1.^2;
            y_square_tavg_acthsf1 = trapz(Rt,y_square_acthsf1,1)./Rt(end);
            y_tavg_acthsf1 = trapz(Rt,Rpop_acthsf1,1)./Rt(end);
            y_tavg_square_acthsf1 = y_tavg_acthsf1.^2;

            ypopavg_acthsf1 = sum(Rpop_acthsf1,2)./length(Rpop_acthsf1(1,:));
            ypopavg_square_acthsf1 = ypopavg_acthsf1.^2;
            ypopavg_square_tavg_acthsf1 = trapz(Rt,ypopavg_square_acthsf1,1)./Rt(end);
            ypopavg_tavg_acthsf1 = trapz(Rt,ypopavg_acthsf1,1)./Rt(end);
            ypopavg_tavg_square_acthsf1 = ypopavg_tavg_acthsf1.^2;

            Rsyn_acthsf1 = (ypopavg_square_tavg_acthsf1 - ypopavg_tavg_square_acthsf1)/(sum((y_square_tavg_acthsf1 - y_tavg_square_acthsf1),2)/length(Rpop_acthsf1(1,:)));

            y_square_mpercry = Rpop_mpercry.^2;
            y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
            y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
            y_tavg_square_mpercry = y_tavg_mpercry.^2;

            ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
            ypopavg_square_mpercry = ypopavg_mpercry.^2;
            ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
            ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
            ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

            Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
            Rsyn = [Rsyn;Rsyn_acthsf1,Rsyn_mpercry];
    end
    
    rn = rn + 1
end

line_color = {[52 57 60]./255,[219 73 76]./255};
figure(2) 
legend(tempamp_name,'Orientation','vertical','Location','southoutside','FontSize',12);

figure %Fig.4C
errorbar([-1,tempamp_test(2:end)],phasedistri_acthsf1(:,1),phasedistri_acthsf1(:,2),'s--','Color',line_color{1},'LineWidth',1.5,'MarkerEdgeColor',line_color{1}); hold on
errorbar([-1,tempamp_test],phasedistri_mpercry(:,1),phasedistri_mpercry(:,2),'o--','Color',line_color{2},'LineWidth',1.5,'MarkerEdgeColor',line_color{2}); hold on
xlabel('Temperature cycle amplitudes (°C)'); ylim([5,20]); ylabel('Phase distribution (mean \pm SD)'); set(gca,'XTick',-1:1:5,'XTickLabel',{tempamp_name{1},'0','1','2','3','4','5'},'XTickLabelRotation',0,'FontSize',16); axis square
legend('actHSF1','Per/Cry mRNA','Location','southeast');

figure %Fig.4D
errorbar([-1,tempamp_test(2:end)],perioddistri_acthsf1(:,1),perioddistri_acthsf1(:,2),'s--','Color',line_color{1},'LineWidth',1.5,'MarkerEdgeColor',line_color{1}); hold on
errorbar([-1,tempamp_test],perioddistri_mpercry(:,1),perioddistri_mpercry(:,2),'o--','Color',line_color{2},'LineWidth',1.5,'MarkerEdgeColor',line_color{2}); hold on
xlabel('Temperature cycle amplitudes (°C)'); ylim([23.6,24.6]); ylabel('Period distribution (mean \pm SD)'); set(gca,'XTick',-1:1:5,'XTickLabel',{tempamp_name{1},'0','1','2','3','4','5'},'XTickLabelRotation',0,'FontSize',16); axis square
legend('actHSF1','Per/Cry mRNA','Location','southeast');

figure %Fig.4E
plot([-1,tempamp_test(2:end)],Rsyn([1,3:end],1),'s--','Color',line_color{1},'LineWidth',1.5); hold on
plot([-1,tempamp_test],Rsyn(:,2),'o--','Color',line_color{2},'LineWidth',1.5); hold on
xlabel('Temperature cycle amplitudes (°C)'); ylim([0,1]); ylabel('R_{syn}'); axis square
legend('actHSF1','Per/Cry mRNA','Location','southeast');
set(gca,'XTick',-1:1:5,'XTickLabel',{tempamp_name{1},'0','1','2','3','4','5'},'XTickLabelRotation',0,'FontSize',16);


%% Supp Fig.S1: average-varying temperature rhythms
clear all; close all; clc;

pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop = pop0;

tRange = 0:1e-1:2400+24*50; %[0 2400+24*50];
Y0 = ones(1*size(pop,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%human body temperature fluctuation
pt = 24;
pw = 12;
Tempint = 3; %°C
p = 0;

tempavg_test = 32:2:40;
tempavg_name = {'T_{avg} = 32 °C','T_{avg} = 34 °C','T_{avg} = 36 °C','T_{avg} = 38 °C','T_{avg} = 40 °C'};
line_color = {[52 57 60]./255, [219 73 76]./255, [235 164 50]./255, [134 170 109]./255, [70 181 211]./255, [102 103 171]./255};

rn = 1
for i = 1:length(tempavg_test)
    
    Tcoreavg = tempavg_test(i);
    tic;
    [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync(t,Y,pw,pt,Tempint,Tcoreavg,pop,p),tRange,Y0,options);
    toc
    YY = reshape(Y,[],11,size(pop,1));
    mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
    mpercry_ens = mean(mpercry_pop,2);
    acthsf1_pop = reshape(YY(:,8,:),[],size(pop,1));
    acthsf1_ens = mean(acthsf1_pop,2);
    
    figure(1) %Fig.S1A
    plot(t-2400,acthsf1_ens,'LineWidth',2,'Color',line_color{i}); hold on; 
    xlim([0 48]); xlabel('Time (h)'); ylabel({('actHSF1'),('ensemble average level (a.u.)')}); 
    set(gca,'XTick',0:6:48,'FontSize',16); box on; axis square
    figure(2) %Fig.S1B
    plot(t-2400,mpercry_ens,'LineWidth',2,'Color',line_color{i}); hold on; 
    xlim([0 48]); xlabel('Time (h)'); ylabel({('Per/Cry mRNA'),('ensemble average level (a.u.)')}); 
    set(gca,'XTick',0:6:48,'FontSize',16); box on; axis square
    
    rn = rn + 1
end
figure(2) 
legend(tempavg_name,'Orientation','vertical','Location','eastoutside','FontSize',14);

%% Fig.S2: Arnold Tongue
%% Computational results
%"brute force integration method" (Schmal et al., 2015)
%clear all; clc; close all;

pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop = pop0;


options = odeset('AbsTol',1e-10,'RelTol',1e-10);
Y0 = ones(1*size(pop,1),11);

%characteristics of temperature rhythm
Tempint_var = 1:1:5; %0:1:5; %entrainer strength/amplitude (°C)
pt_var = 22:0.2:26; %entrainer period (h)
Tcoreavg = 37; %°C
p = 0;

arnoldtongue_Rsyn = [];
arnoldtongue = [];
rn = 1
for i = 1:length(Tempint_var) 
    Tempint = Tempint_var(i);
    rn1 = 1
    for j = 1:length(pt_var)
        pt = pt_var(j);
        pw = 0.5*pt; %keep the thermoperiod same at 50%
        tspan = 0:1.5e-1:pt*100+pt*100; %[0 pt*200];
        
        [t,y] = ode45(@(t,Y)Temp_HSF1_PCG_sync(t,Y,pw,pt,Tempint,Tcoreavg,pop,p),tspan,Y0,options);
        Y = reshape(y,[],11,size(pop,1));
        mpercry_pop = reshape(Y(:,1,:),[],size(pop,1));
        mpercry_ens = mean(mpercry_pop,2);
        
        index0 = find(t>pt*100+pt*84 & t<=pt*100+pt*85); %select the 85th entrainment cycle, 84T < t <= 85T (T is the entrainer period), to determine the reference state y0
        mpercry_ens0 = mpercry_ens(index0);
        indmax0 = find(mpercry_ens0 == max(mpercry_ens0)); %choose the state occurring once per cycle, such as max or min, and of a specific species to ensure the state to be unique per cycle
        y0 = mpercry_ens0(indmax0); 
        %determine the reccurence times t_recc during t > 85T by looking for all y1's that are the same as y0 using the criterion that ||y1 − y0||_2 < ϵ} with ϵ=0.02, ||A||_2 is 2-norm
        index1 = find(t>pt*100+pt*85);
        t1 = t(index1);
        mpercry_ens1 = mpercry_ens(index1); 
        y1_norm = [];
        for n = 1:length(mpercry_ens1)
            n
            y1_norm = [y1_norm,norm((mpercry_ens1(n)-y0),2)];
        end
        index2 = find(diff(sign(diff(y1_norm)))>0)+1;
        y1_norm_min = y1_norm(index2);
        eps = 0.02;
        index3 = find(y1_norm_min<eps);
        t_recc = t1(index2(index3));
        t_recc_remov = []; %eliminate the (systematic) errors caused by software numerical accuracy
        for m = 2:length(t_recc)
            if t_recc(m)-t_recc(m-1) < 5
                t_recc_remov(end+1) = t_recc(m);
            end
        end
        t_recc = setdiff(t_recc,t_recc_remov);
        %determine whether the oscillator (ensemble) gets entrained with two criteria: 
        %(1) the recurrence times do not change overtime (|tn+1 − tn| < δ with δ being small)
        %(2) their differences approach the entrainer period T
        if mean(abs(diff(diff(t_recc)))) < 0.3 & abs(mean(diff(t_recc))-pt) < 0.1
            
            indmax_mpercry_ens = find(diff(sign(diff(mpercry_ens1)))<0)+1;
            pkt_mpercry_ens = t1(indmax_mpercry_ens(1))-(pt*100+pt*85);
            pd_mpercry_ens = mean(diff(t1(indmax_mpercry_ens)));
            pkt_pd_mpercry_ens = mod(pkt_mpercry_ens,pd_mpercry_ens);
             
            arnoldtongue = [arnoldtongue;Tempint,pt,pkt_pd_mpercry_ens];
        end
        
        
        Rpop_mpercry = mpercry_pop(index1,:);
        Rt = t1-(pt*100+pt*85);
        
        y_square_mpercry = Rpop_mpercry.^2;
        y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
        y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
        y_tavg_square_mpercry = y_tavg_mpercry.^2;

        ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
        ypopavg_square_mpercry = ypopavg_mpercry.^2;
        ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
        ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
        ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

        Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
        arnoldtongue_Rsyn = [arnoldtongue_Rsyn;Tempint,pt,Rsyn_mpercry];
        
        save('arnoldtongueResult.mat','arnoldtongue_Rsyn','arnoldtongue');
        rn1 = rn1+1
     end
     rn = rn+1
end
%% Plots
%close all; clear all; clc;

load('arnoldtongueResult.mat')
at_rsyn = arnoldtongue_Rsyn;
at = arnoldtongue;

%Fig.S2A
scatter(at_rsyn(:,2),at_rsyn(:,1),200,at_rsyn(:,3),'filled'); hold on
colorbar('Ticks',0:0.1:1); caxis([0 1]);
xlim([21.8 26.2]); ylim([1 5]); xlabel('Temperature period (h)'); ylabel('Temperature amplitude (°C)'); 
set(gca,'FontSize',16); axis square; box on
title('R_{syn} for Per/Cry mRNA','FontSize',24,'FontName','Times New Roman');

%Fig.S2B,C
index = find(at_rsyn(:,3)<0.4); %not show the results with almost desynchronization
for i = 1:size(at,1)
    for j = 1:length(index)
        if at(i,1) == at_rsyn(index(j),1) & at(i,2) == at_rsyn(index(j),2)
            at(i,:) = nan;
        end
    end
end

figure
scatter(at(:,2),at(:,1),200,at(:,3),'filled'); hold on
colorbar; caxis([10 16]); 
xlim([21.8 26.2]); ylim([1 5]); xlabel('Temperature period (h)'); ylabel('Temperature amplitude (°C)'); 
xticks(22:1:26); set(gca,'FontSize',16); axis square; box on
title({('Per/Cry mRNA ensemble average'),('phase (h)')},'FontSize',24,'FontName','Times New Roman')

figure
scatter(at(:,2),at(:,1),200,at(:,3)-at(:,2)/2,'filled'); hold on
colormap(hot); colorbar;
xlim([21.8 26.2]); ylim([1 5]); xlabel('Temperature period (h)'); ylabel('Temperature amplitude (°C)'); 
xticks(22:1:26); set(gca,'FontSize',16); axis square; box on
title({('Per/Cry mRNA ensemble average'),('phase - Onset of the cold period (h)')},'FontSize',24,'FontName','Times New Roman')

%% Fig.5: alternating shift schedules
clear all; close all; clc;

pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop = pop0;

tRange = 0:2e-1:2400+24*150; %[0 2400+24*150]; 
Y0 = ones(1*size(pop,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%characteristics of temperature rhythm
pt = 24;
pw = 12;
Tcoreavg = 37; %°C
Tempint = 3; %°C

p = 2; %alternating "shiftwork" type schedule
Tempint_sw = Tempint;
swsched_test = [5 7; %5-day inverse per week
                3 7; %3-day inverse per week
                20 28; %20-day inverse per 4-week: 5:2
                12 28]; %12-day shift per 4-week: 3:4               

ttl = {'7-5:2 (R:N) ASS','7-3:4 (R:N) ASS',...
       '28-5:2 (R:N) ASS','28-3:4 (R:N) ASS'};                
flname = {'ode45Output-5-2SW.mat','ode45Output-3-4SW.mat','ode45Output-20-8SW.mat','ode45Output-12-16SW.mat'};
              
Rsyn = cell(size(swsched_test,1),2);
rn = 1
for i = 1:size(swsched_test,1)
    swsched = swsched_test(i,:);
    tic;
    [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync_shiftw(t,Y,pw,pt,Tempint,Tcoreavg,pop,p,swsched,Tempint_sw),tRange,Y0,options);
    toc
    YY = reshape(Y,[],11,size(pop,1));
    mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
    mpercry_ens = mean(mpercry_pop,2);
    
    figure(i); 
    tl = tiledlayout(1,2,'TileSpacing','compact');
    bkgAx = axes(tl,'XTick',[],'YTick',[],'Box','on');
    bkgAx.Layout.TileSpan = [1 2];
    ax1 = axes(tl); 
    pkt_mpercry_pop = []; pkt_mpercry_ensavg = [];
    pn = 1
    for n = 1:size(pop,1)  
        indmax_mpercry = find(diff(sign(diff(mpercry_pop(:,n))))<0)+1;
        pkt_mpercry = t(indmax_mpercry);
        pkt_mpercry_24 = mod(pkt_mpercry,24);
        
        plot(ax1,pkt_mpercry_24, pkt_mpercry-1200,'o','color',[0.5 0.5 0.5],'MarkerSize',2); hold on
        pn = pn + 1
    end
    indmax_mpercry_ens = find(diff(sign(diff(mpercry_ens)))<0)+1;
    pkt_mpercry_ens = t(indmax_mpercry_ens);
    pkt_mpercry_ens_24 = mod(pkt_mpercry_ens,24);

    plot(ax1,pkt_mpercry_ens_24, pkt_mpercry_ens-1200,'ro','MarkerFaceColor','r','MarkerSize',4); hold on; 
    xlim(ax1,[0 24]); ylim(ax1,[0 tRange(end)-1200]); xlabel(ax1,{('Per/Cry mRNA'),('ensemble peaking time (h)')}); ylabel(ax1,'Time (day)'); 
    set(ax1,'XTick',0:6:24,'YDir','reverse','YTick',0:24*50:tRange(end)-1200,'YTickLabel',{'0','50','100','150','200'},'FontSize',16); box on;
            
    
    Rsyn{i,1} = swsched; 
    Rsyn0 = [];
    for interv = 1:t(end)/24
        interv
        
        index = find(t>=24*(interv-1) & t<24*interv);
        Rpop_mpercry = mpercry_pop(index,:);
        Rt = t(index)-24*(interv-1);
        
        y_square_mpercry = Rpop_mpercry.^2;
        y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
        y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
        y_tavg_square_mpercry = y_tavg_mpercry.^2;

        ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
        ypopavg_square_mpercry = ypopavg_mpercry.^2;
        ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
        ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
        ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

        Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
        Rsyn0 = [Rsyn0;Rsyn_mpercry];
     end
     Rsyn{i,2} = Rsyn0;
    
     figure(i); 
     ax2 = axes(tl); ax2.Layout.Tile = 2;
     plot(ax2,Rsyn0,(t(t~=0 & mod(t,24)==0)-1200)/24,'ks','MarkerFaceColor','k','MarkerSize',4); hold on; 
     xlim(ax2,[0 1]); ylim(ax2,[0 (tRange(end)-1200)/24]); xlabel(ax2,'R_{syn} for Per/Cry mRNA'); 
     set(ax2,'XTick',0:0.2:1,'YDir','reverse','YTick',0:24*50:tRange(end)-1200,'YTickLabel',{' ','50','100','150','200','250','300'},'FontSize',16); 
     ax2.XAxisLocation = 'bottom'; ax2.XColor = 'k'; ax2.Color = 'none'; box on; %ax2.YAxis.Visible = 'off';  
     title(tl,ttl{i},'FontName','Times New Roman','FontSize',24);
     
     save(flname{i},'t','Y','Rsyn');
     rn = rn + 1
end
     
%% Fig.6: amplitude-varying alternating shift schedules
clear all; close all; clc;

pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop = pop0;

tRange = 0:2e-1:2400+24*150; %[0 2400+24*150]; 
Y0 = ones(1*size(pop,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%characteristics of temperature rhythm
pt = 24;
pw = 12;
Tcoreavg = 37; %°C
Tempint = 3; %°C

p = 2; %alternating "shiftwork" type schedule
Tempint_sw_test = [5,3,1];
swsched_test = [5 7; %5-day inverse per week
                3 7]; %3-day inverse per week
                
ttl = {'7-5:2 (R(with \DeltaT varied):N) ASS','7-3:4 (R(with \DeltaT varied):N) ASS'};                
flname = {'ode45Output-5-2SW-increasedTamp.mat','ode45Output-5-2SW-reducedTamp.mat',...
          'ode45Output-3-4SW-increasedTamp.mat','ode45Output-3-4SW-reducedTamp.mat'};
lgd = {'\DeltaT_r = 5 °C','\DeltaT_r = 3 °C','\DeltaT_r = 1 °C'};
markertyp = {'rs','ko','bv'}; markercolor = {'r','k','b'};
Rsyn = cell(2*size(swsched_test,1),3);
rn = 1
for i = 1:size(swsched_test,1)
    swsched = swsched_test(i,:);
    rn1 = 1
    figure(i); 
    tl = tiledlayout(1,2,'TileSpacing','compact');
    bkgAx = axes(tl,'XTick',[],'YTick',[],'Box','on');
    bkgAx.Layout.TileSpan = [1 2];
    for j = 1:length(Tempint_sw_test)
        Tempint_sw = Tempint_sw_test(j);
        tic;
        [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync_shiftw(t,Y,pw,pt,Tempint,Tcoreavg,pop,p,swsched,Tempint_sw),tRange,Y0,options);
        toc
        YY = reshape(Y,[],11,size(pop,1));
        mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
        mpercry_ens = mean(mpercry_pop,2);
        indmax_mpercry_ens = find(diff(sign(diff(mpercry_ens)))<0)+1;
        pkt_mpercry_ens = t(indmax_mpercry_ens);
        pkt_mpercry_ens_24 = mod(pkt_mpercry_ens,24);

        ax1 = axes(tl); ax1.Layout.Tile = 1;
        pl1(j)=plot(ax1,pkt_mpercry_ens_24, pkt_mpercry_ens-1200,markertyp{j},'MarkerFaceColor',markercolor{j},'MarkerSize',5); hold on; 
        xlim(ax1,[0 24]); ylim(ax1,[0 tRange(end)-1200]); xlabel(ax1,{('Per/Cry mRNA'),('ensemble peaking time (h)')}); ylabel(ax1,'Time (day)'); 
        ax1.Color = 'none'; box on;
        set(ax1,'XTick',0:6:24,'YDir','reverse','YTick',0:24*50:tRange(end)-1200,'YTickLabel',{'0','50','100','150','200'},'FontSize',16);
        if j==3
            legend(pl1,lgd,'Location','northwest','Color','none','EdgeColor','none','FontName','Times New Roman','FontSize',13);
        end


        Rsyn{i,1} = swsched; Rsyn{i,2} = Tempint_sw;
        Rsyn0 = [];
        for interv = 1:t(end)/24
            interv

            index = find(t>=24*(interv-1) & t<24*interv);
            Rpop_mpercry = mpercry_pop(index,:);
            Rt = t(index)-24*(interv-1);

            y_square_mpercry = Rpop_mpercry.^2;
            y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
            y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
            y_tavg_square_mpercry = y_tavg_mpercry.^2;

            ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
            ypopavg_square_mpercry = ypopavg_mpercry.^2;
            ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
            ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
            ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

            Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
            Rsyn0 = [Rsyn0;Rsyn_mpercry];
         end
         Rsyn{i,3} = Rsyn0;

         ax2 = axes(tl); ax2.Layout.Tile = 2;
         pl2(j)=plot(ax2,Rsyn0,(t(t~=0 & mod(t,24)==0)-1200)/24,markertyp{j},'MarkerFaceColor',markercolor{j},'MarkerSize',5); hold on; 
         xlim(ax2,[0 1]); ylim(ax2,[0 (tRange(end)-1200)/24]); xlabel(ax2,'R_{syn} for Per/Cry mRNA'); 
         set(ax2,'XTick',0:0.2:1,'YDir','reverse','YTick',0:24*50:tRange(end)-1200,'YTickLabel',{' ','50','100','150','200','250','300'},'FontSize',16); 
         ax2.XAxisLocation = 'bottom'; ax2.XColor = 'k'; ax2.Color = 'none'; box on; 
         title(tl,ttl{i},'FontName','Times New Roman','FontSize',24); %,'FontWeight','bold');
         if j==3
            legend(pl2,lgd,'Location','northwest','Color','none','EdgeColor','none','FontName','Times New Roman','FontSize',13);
         end
         
         if i == 1
             if j == 1
                 save(flname{1},'t','Y','Rsyn');
             elseif j == 3
                 save(flname{2},'t','Y','Rsyn');
             end
         elseif i == 2
             if j == 1
                 save(flname{3},'t','Y','Rsyn');
             elseif j == 3
                 save(flname{4},'t','Y','Rsyn');
             end
         end

         rn1 = rn1 + 1
    end
    rn = rn + 1
end

%% Select virtual individual population and Fig.S3A
clear all; close all; clc;

%nominal individual
tempconvt = [0.13, 0.65, 0.65,... %kin0, kin1, kout1
             8, 20, 0.45]; %kin2, Kin2, kout2 
%construct individual population         
npop=300;
s = sobolset(length(tempconvt));
varmin = [0.05, 0.4, 0.65, 8, 5,  0.45]; %different combinations of values of kin0, kin1, Kin2 represent individual differences
varmax = [0.5,  1,   0.65, 8, 30, 0.45];
for n = 1:npop
    pop00 = varmin+(varmax-varmin).*s(n,:);
    pop0(n,:) = pop00;
end
writematrix(pop0,'pop_TS_samp.xlsx','WriteMode','append');

%select valid individuals 
pop1 = xlsread('pop_ALL_samp3%.xlsx');
pop1 = pop1(500:end,7:end); %assume each individual has 500 cells, which are the same in terms of HSF1-involved HSR pathway and core clock for individuals
  
tRange = 0:1e-1:1200+24*70; %[0 1200+24*70];
Y0 = ones(1*size(pop1,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

pt = 24;
pw = 12;
Tcoreavg = 37; %°C
Tempint = 3; %°C
p = 0;
for i = 1:size(pop0,1)
    i
    indvd0 = [repmat(pop0(i,:),size(pop1,1),1),pop1];
    pop = indvd0;
    tic;
    [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync(t,Y,pw,pt,Tempint,Tcoreavg,pop,p),tRange,Y0,options);
    toc
    YY = reshape(Y,[],11,size(pop,1));
    mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
    mpercry_ens = mean(mpercry_pop,2);
    
    index = find(t>=2400);
    t_period = t(index);
    indmax_mpercry_ens = find(diff(sign(diff(mpercry_ens(index))))<0)+1;
    period_mpercry_ens = mean(diff(t_period(indmax_mpercry_ens)));

    Rpop_mpercry = mpercry_pop(index,:);
    Rt = t(index)-2400;
    
    y_square_mpercry = Rpop_mpercry.^2;
    y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
    y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
    y_tavg_square_mpercry = y_tavg_mpercry.^2;

    ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
    ypopavg_square_mpercry = ypopavg_mpercry.^2;
    ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
    ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
    ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

    Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));

    %selection criteria: under nominal/normal T rhythm, the host (1) gets entrained with a 24-h period, and (2) has an almost full synchronization of its cell population
    if abs(period_mpercry_ens-24)<0.02 & Rsyn_mpercry>=0.8
        scatter3(pop0(i,1),pop0(i,2),pop0(i,5),35,Rsyn_mpercry,'filled'); hold on; 
        xlabel('$v_{act0,ts1}$','Interpreter','latex'); ylabel('$v_{act1,ts1}$','Interpreter','latex'); zlabel('$k_{b,ts2}$','Interpreter','latex'); 
        h = colorbar; h.Ticks = 0.84:0.02:0.98; h.FontSize = 16;
        title('Mean R_{syn,nom}','FontName','Times New Roman','FontWeight','bold'); set(gca,'FontSize',16); box on; axis square;
        
        pop0_select(i,:) = [pop0(i,:),Rsyn_mpercry];
        writematrix(pop0_select(i,:),'pop0_select','WriteMode','append')
    end
end

%% Fig.7 and S3B: individual internal synchrony degree under ASSs
%% Computational results
close all; clc; clear all; 

pop_select = load('pop0_select.txt');
pop1 = pop_select(:,1:6);
pop0 = xlsread('pop_ALL_samp3%.xlsx');
pop2 = pop0(1:500,7:end); %each individual has 500 cells, which are the same in terms of HSF1-involved HSR pathway and core clock for individuals

tRange = 0:1.5e-1:2400+24*150; %[0 2400+24*150];
Y0 = ones(1*size(pop2,1),11);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%characteristics of temperature alternating shift schedules
pt = 24;
pw = 12;
Tcoreavg = 37; %°C
Tempint = 3; %°C
p = 2;
Tempint_sw = Tempint;
swsched_sw = [6 7; 5 7; 3 7];
flname = {'shiftworkexpResult6-1.mat','shiftworkexpResult5-2.mat','shiftworkexpResult3-4.mat'};

sw_indivdiff = cell(size(pop1,1),3); %1st col, sw schedule; 2nd col, ensemble average peaking time; 3rd col, Rsyn
for n = 1:size(swsched_sw,1)
    n
    swsched = swsched_sw(n,:);
    for i = 1:size(pop1,1)
        i
        indvd0 = [repmat(pop1(i,:),size(pop2,1),1),pop2];
        pop = indvd0;
        tic;
        [t,Y] = ode45(@(t,Y)Temp_HSF1_PCG_sync_shiftw(t,Y,pw,pt,Tempint,Tcoreavg,pop,p,swsched,Tempint_sw),tRange,Y0,options);
        toc
        YY = reshape(Y,[],11,size(pop,1));
        mpercry_pop = reshape(YY(:,1,:),[],size(pop,1));
        mpercry_ens = mean(mpercry_pop,2);

        indmax_mpercry_ens = find(diff(sign(diff(mpercry_ens)))<0)+1;
        pkt_mpercry_ens = t(indmax_mpercry_ens);
        pkt_mpercry_ens_24 = mod(pkt_mpercry_ens,24);
        
        sw_indivdiff{i,1} = swsched;
        sw_indivdiff{i,2} = pkt_mpercry_ens;

        Rsyn0 = [];
        for interv = 1:t(end)/24
            interv
            index = find(t>=24*(interv-1) & t<24*interv);

            Rpop_mpercry = mpercry_pop(index,:);
            Rt = t(index)-24*(interv-1);
            
            y_square_mpercry = Rpop_mpercry.^2;
            y_square_tavg_mpercry = trapz(Rt,y_square_mpercry,1)./Rt(end);
            y_tavg_mpercry = trapz(Rt,Rpop_mpercry,1)./Rt(end);
            y_tavg_square_mpercry = y_tavg_mpercry.^2;

            ypopavg_mpercry = sum(Rpop_mpercry,2)./length(Rpop_mpercry(1,:));
            ypopavg_square_mpercry = ypopavg_mpercry.^2;
            ypopavg_square_tavg_mpercry = trapz(Rt,ypopavg_square_mpercry,1)./Rt(end);
            ypopavg_tavg_mpercry = trapz(Rt,ypopavg_mpercry,1)./Rt(end);
            ypopavg_tavg_square_mpercry = ypopavg_tavg_mpercry.^2;

            Rsyn_mpercry = (ypopavg_square_tavg_mpercry - ypopavg_tavg_square_mpercry)/(sum((y_square_tavg_mpercry - y_tavg_square_mpercry),2)/length(Rpop_mpercry(1,:)));
            Rsyn0 = [Rsyn0,Rsyn_mpercry];
        end
            sw_indivdiff{i,3} = Rsyn0;
    end
    save(flname{n},'sw_indivdiff')
end
%% Plots
clc; clear all; close all; 

pop_select = load('pop0_select.txt');
pop1 = pop_select(:,1:6);

swdata = cell(size(pop1,1),3); 
statdata = cell(3,3);
markclr = {'r','k','b'};
lgd = {'7-6:1 (R:N) ASS','7-5:2 (R:N) ASS','7-3:4 (R:N) ASS'};
marksize = 60;
for n = 1:3
    if n==1
        load('shiftworkexpResult6-1.mat');
        swdata = sw_indivdiff;
    elseif n==2
        load('shiftworkexpResult5-2.mat');
        swdata = sw_indivdiff;
    else
        load('shiftworkexpResult3-4.mat');
        swdata = sw_indivdiff;
    end
    
    tcal_sw = 100;
    trange = 0:1.5e-1:2400+24*150;
    t = trange(trange~=0 & mod(trange,24)==0); %24/0.15=160
    index_nom = find(t>1200 & t<=1200+24*50);
    index_sw = find(t>1200+24*tcal_sw);
    sw_indivdiff_index_rsyn = []; %rsynmean_nom, rsyndelta_nom, rsynmean_sw, rsyndelta_sw
    for i = 1:size(pop1,1)
        i
        rsyn = cell2mat(swdata(i,3));
        rsyn_nom = rsyn(index_nom);
        rsyn_sw = rsyn(index_sw);
        
        rsynmean_nom = mean(rsyn_nom);
        rsyndelta_nom = max(rsyn_nom)-min(rsyn_nom);
        rsynmean_sw = mean(rsyn_sw);
        rsyndelta_sw = max(rsyn_sw)-max(rsyn_sw);
        sw_indivdiff_index_rsyn = [sw_indivdiff_index_rsyn; rsynmean_nom,rsyndelta_nom,rsynmean_sw,rsyndelta_sw]; 
    end
    
    %Statistical analysis (linear regression) and Fig.S3B
    statdata{n,1} = swdata{1,1};
    mdl_nom = fitlm([pop1(:,1),pop1(:,2),pop1(:,5)],sw_indivdiff_index_rsyn(:,1));
    mdl_sw = fitlm([pop1(:,1),pop1(:,2),pop1(:,5)],sw_indivdiff_index_rsyn(:,3));
    statdata{n,2} = mdl_nom;
    statdata{n,3} = mdl_sw;
   
    if n~=1
        figure(1)
        sgtitle('Mean R_{syn,ASS}','FontSize',20,'FontWeight','bold','FontName','Times New Roman');
        subplot(1,2,n-1)
        scatter3(pop1(:,1),pop1(:,2),pop1(:,5),marksize,sw_indivdiff_index_rsyn(:,3),'filled'); hold on; 
        h = colorbar; h.FontSize = 16;
        if n==2
            h.Ticks = 0.3:0.1:0.9; 
        else
            h.Ticks = 0.1:0.1:0.8;
        end
        xlabel('$v_{act0,ts1}$','Interpreter','latex'); ylabel('$v_{act1,ts1}$','Interpreter','latex'); zlabel('$k_{b,ts2}$','Interpreter','latex'); 
        title(lgd{n},'FontName','Times New Roman'); set(gca,'FontSize',16); axis square; box on; 
    end
        
    %Curve fit and Fig.7
    if n~=3
        f = fit(sw_indivdiff_index_rsyn(:,1),sw_indivdiff_index_rsyn(:,3),'rat15','StartPoint',-ones(1,7));
    else
        ft = fittype( 'p1+p2/(1+exp(-(x-p3)/p4))^p5','coefficients',{'p1','p2','p3','p4','p5'});
        bnd = [0.04,1,0.97,0.012,1];
        options1 = fitoptions('Method','NonlinearLeastSquares','Lower',bnd,'Upper',bnd,'Startpoint',ones(1,5));
        f = fit(sw_indivdiff_index_rsyn(:,1),sw_indivdiff_index_rsyn(:,3),ft,options1);
    end

    figure(2)
    pl(n)=scatter(sw_indivdiff_index_rsyn(:,1),sw_indivdiff_index_rsyn(:,3),marksize,'LineWidth',1.5,'MarkerEdgeColor',markclr{n},'MarkerFaceColor',markclr{n}); alpha 0.45; hold on;
    x_cf = 0.8:(1-0.8)/100:1;
    plot(x_cf,f(x_cf),'LineWidth',1,'Color',markclr{n},'LineStyle','--'); hold on;
    xlim([0.8 1]); ylim([0 1]); xlabel('Mean R_{syn,nom}'); ylabel('Mean R_{syn,ASS}');
    set(gca,'FontSize',20); axis square; box on; grid on; 
end
figure(2); legend(pl,lgd,'Location','eastoutside','Orientation','vertical','FontSize',18,'FontName','Times New Roman');

%% ODE function_single cell
function dYdt = Temp_HSF1_PCG(t,Y,pw,pt,Tempint,Tcoreavg,pcg,hs,entrainparam)
         
         PerCry_mRNA = Y(1); PerCry = Y(2); PerCry_nuc = Y(3); 
         Bmal1_mRNA = Y(4); Bmal1 = Y(5); Bmal1_nuc = Y(6);
         ClockBmal1 = Y(7);
         actHSF1 = Y(8); indHSP = Y(9);
         TS1 = Y(10); TS2 = Y(11);
         
         Temp = square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg; %nominal/normal rhythm
         
         %Temperature input and conversion
         kin0 = 0.13; kin1 = 0.65; kout1 = 0.65;
         kin2 = 8; Kin2 = 20; n = 3; kout2 = 0.45;
         
         dTS1dt = kin0*Tcoreavg + kin1*(Temp-Tcoreavg) - kout1*TS1;
         dTS2dt = kin2*TS1^n/(Kin2^n+TS1^n)-kout2*TS2;
         
         ks = entrainparam(1); Ks = entrainparam(2); 
         khsf = entrainparam(3); Khsf = entrainparam(4);
         
         %HSF1-involved HSR dynamics
         HSF1tot = hs(1); vact = hs(2); vina = hs(3);
         vbhsp = hs(4); vdhsp = hs(5);
         
         dactHSF1dt = vact*(HSF1tot-actHSF1)*(1+ks*TS2/(Ks+TS2))/indHSP-vina*actHSF1;
         dindHSPdt = vbhsp*actHSF1/indHSP-vdhsp*indHSP;

         %PCG dynamics and parameter values (Mavroudis et al, 2014)
         v1b = pcg(1); c = 0.01; k1b = pcg(2); k1i = pcg(3); p = 8; k1d = pcg(4); kc = 0; DRN = 1;
         k2b = pcg(5); q = 2; k2d = pcg(6); k2t = pcg(7); 
         k3t = pcg(8); k3d = pcg(9);
         v4b = pcg(10); k4b = pcg(11); r = 3; k4d = pcg(12);
         k5b = pcg(13); k5d = pcg(14); k5t = pcg(15); 
         k6t = pcg(16); k6d = pcg(17); k6a = pcg(18);
         k7a = pcg(19); k7d = pcg(20);

         dPerCry_mRNAdt = v1b*(ClockBmal1+c)*(1+khsf*actHSF1/(Khsf+actHSF1+indHSP/ClockBmal1))/(k1b*(1+(PerCry_nuc/k1i)^p+ClockBmal1+c))-k1d*PerCry_mRNA+kc*DRN/ClockBmal1;
         dPerCrydt = k2b*PerCry_mRNA^q-k2d*PerCry-k2t*PerCry+k3t*PerCry_nuc;
         dPerCry_nucdt = k2t*PerCry-k3t*PerCry_nuc-k3d*PerCry_nuc;
         dBmal1_mRNAdt = v4b*PerCry_nuc^r/(k4b^r+PerCry_nuc^r)-k4d*Bmal1_mRNA;
         dBmal1dt = k5b*Bmal1_mRNA-k5d*Bmal1-k5t*Bmal1+k6t*Bmal1_nuc;
         dBmal1_nucdt = k5t*Bmal1-k6t*Bmal1_nuc-k6d*Bmal1_nuc-k6a*Bmal1_nuc+k7a*ClockBmal1;
         dClockBmal1dt = k6a*Bmal1_nuc-k7a*ClockBmal1-k7d*ClockBmal1;

         dYdt = [dPerCry_mRNAdt;dPerCrydt;dPerCry_nucdt;dBmal1_mRNAdt;dBmal1dt;dBmal1_nucdt;dClockBmal1dt;dactHSF1dt;dindHSPdt;dTS1dt;dTS2dt];
end

%% ODE function_cell population synchronization
function dYdt = Temp_HSF1_PCG_sync(t,Y,pw,pt,Tempint,Tcoreavg,pop,p)

        if p == 0 %nominal/normal
            Temp = square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg;
        elseif p == 1 %constant + nominal/normal + constant
            Temp = Tcoreavg.*(t<=1200+24*10)+(square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg).*(t>1200+24*10 & t<=1200+24*65)+Tcoreavg.*(t>1200+24*65);
        elseif p == -1 % constant + inversed +constant
            Temp = Tcoreavg.*(t<=1200+24*10)+(-square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg).*(t>1200+24*10 & t<=1200+24*65)+Tcoreavg.*(t>1200+24*65);    
        end
            
        for indvd = 1:size(pop,1)
             PerCry_mRNA = Y((11*indvd-11)+1,:); PerCry = Y((11*indvd-11)+2,:); PerCry_nuc = Y((11*indvd-11)+3,:); 
             Bmal1_mRNA = Y((11*indvd-11)+4,:); Bmal1 = Y((11*indvd-11)+5,:); Bmal1_nuc = Y((11*indvd-11)+6,:);
             ClockBmal1 = Y((11*indvd-11)+7,:);
             actHSF1 = Y((11*indvd-11)+8,:); indHSP = Y((11*indvd-11)+9,:);
             TS1 = Y((11*indvd-11)+10,:); TS2 = Y((11*indvd-11)+11,:);
             
             
             tempconvt = pop(indvd,1:6);
             kin0 = tempconvt(1); kin1 = tempconvt(2); kout1 = tempconvt(3); 
             kin2 = tempconvt(4); Kin2 = tempconvt(5); n = 3; kout2 = tempconvt(6); 

             entrainparam = pop(indvd,7:10);
             ks = entrainparam(1); Ks = entrainparam(2); 
             khsf = entrainparam(3); Khsf = entrainparam(4);
             
             hs = pop(indvd,11:15);
             HSF1tot = hs(1); vact = hs(2); vina = hs(3);
             vbhsp = hs(4); vdhsp = hs(5);
             
             pcg = pop(indvd,16:end);
             v1b = pcg(1); c = 0.01; k1b = pcg(2); k1i = pcg(3); p = 8; k1d = pcg(4); kc = 0; DRN = 1;
             k2b = pcg(5); q = 2; k2d = pcg(6); k2t = pcg(7); 
             k3t = pcg(8); k3d = pcg(9);
             v4b = pcg(10); k4b = pcg(11); r = 3; k4d = pcg(12);
             k5b = pcg(13); k5d = pcg(14); k5t = pcg(15); 
             k6t = pcg(16); k6d = pcg(17); k6a = pcg(18);
             k7a = pcg(19); k7d = pcg(20);
             
             
             dTS1dt(indvd) = kin0*Tcoreavg + kin1*(Temp-Tcoreavg) - kout1*TS1;
             dTS2dt(indvd) = kin2*TS1^n/(Kin2^n+TS1^n)-kout2*TS2;

             dactHSF1dt(indvd) = vact*(HSF1tot-actHSF1)*(1+ks*TS2/(Ks+TS2))/indHSP-vina*actHSF1;
             dindHSPdt(indvd) = vbhsp*actHSF1/indHSP-vdhsp*indHSP;

             dPerCry_mRNAdt(indvd) = v1b*(ClockBmal1+c)*(1+khsf*actHSF1/(Khsf+actHSF1+indHSP/ClockBmal1))/(k1b*(1+(PerCry_nuc/k1i)^p+ClockBmal1+c))-k1d*PerCry_mRNA+kc*DRN/ClockBmal1;
             dPerCrydt(indvd) = k2b*PerCry_mRNA^q-k2d*PerCry-k2t*PerCry+k3t*PerCry_nuc;
             dPerCry_nucdt(indvd) = k2t*PerCry-k3t*PerCry_nuc-k3d*PerCry_nuc;
             dBmal1_mRNAdt(indvd) = v4b*PerCry_nuc^r/(k4b^r+PerCry_nuc^r)-k4d*Bmal1_mRNA;
             dBmal1dt(indvd) = k5b*Bmal1_mRNA-k5d*Bmal1-k5t*Bmal1+k6t*Bmal1_nuc;
             dBmal1_nucdt(indvd) = k5t*Bmal1-k6t*Bmal1_nuc-k6d*Bmal1_nuc-k6a*Bmal1_nuc+k7a*ClockBmal1;
             dClockBmal1dt(indvd) = k6a*Bmal1_nuc-k7a*ClockBmal1-k7d*ClockBmal1;
          
        end
        dYdt = [dPerCry_mRNAdt;dPerCrydt;dPerCry_nucdt;dBmal1_mRNAdt;dBmal1dt;dBmal1_nucdt;dClockBmal1dt;dactHSF1dt;dindHSPdt;dTS1dt;dTS2dt];
        dYdt = dYdt(:);            
end

%% ODE function_synchronization under "shift work" type temperature schedule
function dYdt = Temp_HSF1_PCG_sync_shiftw(t,Y,pw,pt,Tempint,Tcoreavg,pop,p,swsched,Tempint_sw)

        if p == 2 %alternating "shiftwork": e.g., 5-day inverse per week (swsched = [5,7]), 14-day inverse per 4-week (swsched = [14,28])
            if mod(t/24,swsched(2)) <= swsched(1)
                Temp_swsched = -square(t*(2*pi)/pt,(pw/pt)*100)*Tempint_sw/2+Tcoreavg;
            elseif mod(t/24,swsched(2)) > swsched(1)
                Temp_swsched = square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg;
            end
            Temp = (square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg).*(t<=2400)...
                 + Temp_swsched.*(t>2400);
        end
             
        for indvd = 1:size(pop,1)

             PerCry_mRNA = Y((11*indvd-11)+1,:); PerCry = Y((11*indvd-11)+2,:); PerCry_nuc = Y((11*indvd-11)+3,:); 
             Bmal1_mRNA = Y((11*indvd-11)+4,:); Bmal1 = Y((11*indvd-11)+5,:); Bmal1_nuc = Y((11*indvd-11)+6,:);
             ClockBmal1 = Y((11*indvd-11)+7,:);
             actHSF1 = Y((11*indvd-11)+8,:); indHSP = Y((11*indvd-11)+9,:);
             TS1 = Y((11*indvd-11)+10,:); TS2 = Y((11*indvd-11)+11,:);
                        
 
             tempconvt = pop(indvd,1:6);
             kin0 = tempconvt(1); kin1 = tempconvt(2); kout1 = tempconvt(3); 
             kin2 = tempconvt(4); Kin2 = tempconvt(5); n = 3; kout2 = tempconvt(6);

             entrainparam = pop(indvd,7:10);
             ks = entrainparam(1); Ks = entrainparam(2); 
             khsf = entrainparam(3); Khsf = entrainparam(4);

             hs = pop(indvd,11:15);
             HSF1tot = hs(1); vact = hs(2); vina = hs(3);
             vbhsp = hs(4); vdhsp = hs(5);
             
             pcg = pop(indvd,16:end); 
             v1b = pcg(1); c = 0.01; k1b = pcg(2); k1i = pcg(3); p = 8; k1d = pcg(4); kc = 0; DRN = 1;
             k2b = pcg(5); q = 2; k2d = pcg(6); k2t = pcg(7); 
             k3t = pcg(8); k3d = pcg(9);
             v4b = pcg(10); k4b = pcg(11); r = 3; k4d = pcg(12);
             k5b = pcg(13); k5d = pcg(14); k5t = pcg(15); 
             k6t = pcg(16); k6d = pcg(17); k6a = pcg(18);
             k7a = pcg(19); k7d = pcg(20);
             
  
             dTS1dt(indvd) = kin0*Tcoreavg + kin1*(Temp-Tcoreavg) - kout1*TS1;
             dTS2dt(indvd) = kin2*TS1^n/(Kin2^n+TS1^n)-kout2*TS2;
  
             dactHSF1dt(indvd) = vact*(HSF1tot-actHSF1)*(1+ks*TS2/(Ks+TS2))/indHSP-vina*actHSF1;
             dindHSPdt(indvd) = vbhsp*actHSF1/indHSP-vdhsp*indHSP;

             dPerCry_mRNAdt(indvd) = v1b*(ClockBmal1+c)*(1+khsf*actHSF1/(Khsf+actHSF1+indHSP/ClockBmal1))/(k1b*(1+(PerCry_nuc/k1i)^p+ClockBmal1+c))-k1d*PerCry_mRNA+kc*DRN/ClockBmal1;
             dPerCrydt(indvd) = k2b*PerCry_mRNA^q-k2d*PerCry-k2t*PerCry+k3t*PerCry_nuc;
             dPerCry_nucdt(indvd) = k2t*PerCry-k3t*PerCry_nuc-k3d*PerCry_nuc;
             dBmal1_mRNAdt(indvd) = v4b*PerCry_nuc^r/(k4b^r+PerCry_nuc^r)-k4d*Bmal1_mRNA;
             dBmal1dt(indvd) = k5b*Bmal1_mRNA-k5d*Bmal1-k5t*Bmal1+k6t*Bmal1_nuc;
             dBmal1_nucdt(indvd) = k5t*Bmal1-k6t*Bmal1_nuc-k6d*Bmal1_nuc-k6a*Bmal1_nuc+k7a*ClockBmal1;
             dClockBmal1dt(indvd) = k6a*Bmal1_nuc-k7a*ClockBmal1-k7d*ClockBmal1;
          
        end
        dYdt = [dPerCry_mRNAdt;dPerCrydt;dPerCry_nucdt;dBmal1_mRNAdt;dBmal1dt;dBmal1_nucdt;dClockBmal1dt;dactHSF1dt;dindHSPdt;dTS1dt;dTS2dt];
        dYdt = dYdt(:);
             
end