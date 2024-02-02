%% STEP1: Manual parameter estimation - for parameters associated with temperature signaling pathway (TS1,TS2)
clear all; close all; clc; 

%characteristic of temperature rhythm (nominal/normal)
pw = 12;
pt = 24;
Tempint = 3;
Tcoreavg = 37;

tspan = [0 480];
y0 = [1,1];
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
A = ode45(@(t,y)tempconv(t,y,pw,pt,Tempint,Tcoreavg),tspan,y0,options);

subplot(1,3,1) %T cycle
temp = square(A.x*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg;
plot(A.x-240,temp,'linewidth',1); hold on
xlim([0 48]); xlabel('Time (h)'); title('temperature (°C)'); axis square
set(gca,'FontSize',14)
subplot(1,3,2) %TS1 (y1) oscillation
plot(A.x-240,A.y(1,:),'linewidth',1); hold on
xlim([0 48]); xlabel('Time (h)'); title('TS1 conc (a.u.)'); axis square
set(gca,'FontSize',14)
subplot(1,3,3) %TS2 (y2) oscillation
plot(A.x-240,A.y(2,:),'linewidth',1); hold on
xlim([0 48]); xlabel('Time (h)'); title('TS2 conc (a.u.)'); axis square
set(gca,'FontSize',14)

index = find(A.x>=240);
amp_y1 = max(A.y(1,index))-min(A.y(1,index));
avg_y1 = (max(A.y(1,index))+min(A.y(1,index)))/2;
disp('TS1 amp/TS1 avg'), disp(amp_y1/avg_y1)

amp_y2 = max(A.y(2,index))-min(A.y(2,index));
avg_y2 = (max(A.y(2,index))+min(A.y(2,index)))/2;
disp('TS2 amp/TS2 avg'), disp(amp_y2/avg_y2)

%acceptance criteria: 
%(1) concentrations of TS1 and TS2 remain positive, and 
%(2) temperature signal is amplified when being transduced through TS1 to TS2, quantified by TS2_amp/TS2_avg > TS1_amp/TS1_avg


%% STEP2: Systematic parameter estimation using the optimization routine - for parameters associated with HSF1-involved HSR pathway and HSF1 entraining effect
%function [pval, fval, exitflag, history, x_initial] = ParamEstim_Temp_HSF1_PCG(~)
clear all; close all; %clc; 

global pt pw Tempint Tcoreavg t_exp PER_exp history1 history2 history3

%data used for calibration
data_exp_cosinor = xlsread('data_exp_cosinor.xlsx'); %post-cosinor data of Fig.1C (Saini et al, 2012)
data_exp_cosinor = data_exp_cosinor';
t_exp = data_exp_cosinor(1,:)+4800; 
PER_exp = data_exp_cosinor(2,:); 

%characteristics of temperature rhythm
pt = 24;
pw = 12;
Tempint = 3; %°C
Tcoreavg = 37; %°C

%main optimization function
lb = [1     1      1      1      5     0.01    0.01   0.1    0.01]; 
ub = [40    20     40     20     30    1       30     10      10];   
x0 = [];
for k=1:length(lb)
    x0 = [x0, rand*(ub(k)-lb(k))+lb(k)];
end
x_initial = x0;


IterSA=200;
IterFS=100;
IterFM=10;
history1 = [];
history2 = [];
history3 = [];

objfun = @obj_PramEstim;
objfun_sa = @obj_sa_PramEstim;

options = saoptimset('Display', 'iter', 'PlotFc', @saplotbestf, 'MaxIter', IterSA,'OutputFcn', @outfun1);
[pval,fval] = simulannealbnd(objfun_sa,x0,lb,ub,options);
'Solution from Simulated Annealing', display(pval), display(fval)
x0 = pval;

options = optimset('Display' , 'iter', 'PlotFcns', @optimplotfval, 'maxIter', IterFS,'OutputFcn', @outfun2,'TolCon',1e-8,'TolX',1e-8); %TolX-StepTolerance
[pval,fval,exitflag,output]=fminsearchcon(objfun,x0,lb,ub,[],[],@constrt_PramEstim,options)
'Solution from direct search ', pval, fval
x0 = pval;

'*** Finalize with gradient optimization'
constrcall_count = 0;
options = optimoptions('fmincon','PlotFcn','optimplotfval','Display','iter-detailed', 'MaxIterations',IterFM, ...
          'Algorithm', 'sqp','OutputFcn', @outfun3,'TolCon',1e-8,'StepTolerance',1e-8); 
[pval,fval,exitflag,output]=fmincon(objfun,x0,[],[],[],[],lb,ub,@constrt_PramEstim,options)
disp(pval)

        
%output functions
history = [history1;history2;history3];
function [stop,options,optchanged]=outfun1(options,optimvalues,state)
global history1
    optchanged = false;
    stop = false;
    if isequal(state,'iter')
    history1 = [history1;optimvalues.x];
    end
end
function [stop,options,optchanged]=outfun2(pval,options,state)
global history2
    optchanged = false;
    stop = false;
    if isequal(state,'iter')
    history2 = [history2;pval];
    end
end
function [stop,options,optchanged]=outfun3(pval,options,state)
global history3
    optchanged = false;
    stop = false;
    if isequal(state,'iter')
    history3 = [history3;pval];
    end
end
%end



%multiple dimensional objective function for SA (objective(s) + nonlinear constraint(s))
function z = obj_sa_PramEstim(x)

%objective
global pt pw Tempint Tcoreavg t_exp PER_exp record0

record0 = [];

tRange = [0 4800+24*50];
Y0 = 1*ones(11,1);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
entrainparam = x(1:4);
hsparam = x(5:end);

A = ode45(@(t,Y)Temp_HSF1_PCG(t,Y,pw,pt,Tempint,Tcoreavg,hsparam,entrainparam),tRange,Y0,options);
index0 = find(A.x>=4800 & A.x<4825);
pkt_expPER2 = t_exp(find(PER_exp == max(PER_exp)));
t0 = A.x(index0);
PERCRY = A.y(2,index0);
pkt_predPERCRY = t0(find(PERCRY == max(PERCRY)));
z1 = (pkt_predPERCRY-pkt_expPER2)^2;
record0 = [record0 z1];

index1 = find(A.x>=4806 & A.x<=4830);
acthsf1 = A.y(8,index1);
amp_acthsf1 = max(acthsf1)-min(acthsf1);
avg_acthsf1 = (max(acthsf1)+min(acthsf1))/2;
z2 = (amp_acthsf1/avg_acthsf1-0.6)^2;
z3 = (avg_acthsf1-1)^2;
record0 = [record0 z2 z3];

z = z1+(z2+z3)*50;
record0 = [record0 z];

%nonlinear constraints implemented using penalty method
z = z+getnonlinear(x);
end

function zz = getnonlinear(x)
zz = 0;
[c,ceq] = constrt_PramEstim(x);

%penalty constant (lambda)
lambda = 10^15; lambdaeq = 10^15;
%inequality constraints (when c=[], length->0)
for k = 1:length(c)
    zz = zz+lambda*c(k)^2*getH(c(k));
end
%equality constraints (when ceq=[], length->0)
for k = 1:length(ceq)
   zz = zz+lambdaeq*ceq(k)^2*geteqH(ceq(k));
end

%test whether inequalities hold
function H = getH(c)
if c <= 0
    H=0; 
else
    H=1; 
end
end
%test whether equalities hold
function H = geteqH(ceq)
if ceq == 0
    H=0;
else
    H=1; 
end
end
end



%general objective function
function sse = obj_PramEstim(x) %SSE: sum of squares error

global pt pw Tempint Tcoreavg t_exp PER_exp record0

record0 = [];
record = [];

tRange = [0 4800+24*50];
Y0 = 1*ones(11,1);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
entrainparam = x(1:4);
hsparam = x(5:end);

A = ode45(@(t,Y)Temp_HSF1_PCG(t,Y,pw,pt,Tempint,Tcoreavg,hsparam,entrainparam),tRange,Y0,options);
index0 = find(A.x>=4800 & A.x<4825);
pkt_expPER2 = t_exp(find(PER_exp == max(PER_exp)));
t0 = A.x(index0);
PERCRY = A.y(2,index0);
pkt_predPERCRY = t0(find(PERCRY == max(PERCRY)));
sse1 = (pkt_predPERCRY-pkt_expPER2)^2;
record0 = [record0 sse1];

index1 = find(A.x>=4806 & A.x<=4830);
acthsf1 = A.y(8,index1);
amp_acthsf1 = max(acthsf1)-min(acthsf1);
avg_acthsf1 = (max(acthsf1)+min(acthsf1))/2;
sse2 = (amp_acthsf1/avg_acthsf1-0.6)^2;
sse3 = (avg_acthsf1-1)^2;
record0 = [record0 sse2 sse3];

sse = sse1+(sse2+sse3)*50;
record0 = [record0 sse];

writematrix(record,'record.txt','WriteMode','append');
end



%constraint function
function [c,ceq] = constrt_PramEstim(x)

%inequality constraints, c
global pt pw Tempint Tcoreavg record0

tRange = [0 4800+24*50];
Y0 = 1*ones(11,1);
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
entrainparam = x(1:4);
hsparam = x(5:end);

A = ode45(@(t,Y)Temp_HSF1_PCG(t,Y,pw,pt,Tempint,Tcoreavg,hsparam,entrainparam),tRange,Y0,options);
index2 = find(A.x>=4800);
t2 = A.x(index2);
for i = 1:7
    compnt = A.y(i,index2);
    mesor = (max(compnt) + min(compnt))/2;
    indmax_period0 = [];
    for j = 2:length(compnt)-1
        if compnt(j-1) < compnt(j) & compnt(j+1) < compnt(j) & mesor < compnt(j)
            indmax_period0(end+1) = j;
        end
    end
    %indmax_period0 = find(diff(sign(diff(A.y(i,index2))))<0)+1;
    t_period0 = t2(indmax_period0);
    t_period_avg(i) = mean(diff(t_period0));
end

eps = 1e-3;
c = sum((t_period_avg-24).^2)-eps;

%equality constraints, ceq
ceq = [];


record = [];
if size(record0,2) < 5 
    record0 = [record0 c+eps];
else
    record0(end) = c+eps;
end
record = [record; record0];
writematrix(record,'record.txt','WriteMode','append');
end


%% ODE functions
%STEP1 function 
function dydt = tempconv(t,y,pw,pt,Tempint,Tcoreavg)
    
    y1 = y(1); y2 = y(2); %y1-TS1, y2-TS2

    Temp = square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg;
    
    kin0 = 0.13; kin1 = 0.65; kout1 = 0.65; 
    kin2 = 8; Kin2 = 20; n = 3; kout2 = 0.45;
    
    dy1dt = kin0*Tcoreavg + kin1*(Temp-Tcoreavg) - kout1*y1;
    dy2dt = kin2*y1^n/(Kin2^n+y1^n)-kout2*y2;
    
    dydt = [dy1dt;dy2dt];
end

%STEP2 function 
function dYdt = Temp_HSF1_PCG(t,Y,pw,pt,Tempint,Tcoreavg,hsparam,entrainparam)

         PerCry_mRNA = Y(1); PerCry = Y(2); PerCry_nuc = Y(3); 
         Bmal1_mRNA = Y(4); Bmal1 = Y(5); Bmal1_nuc = Y(6);
         ClockBmal1 = Y(7);
         actHSF1 = Y(8); indHSP = Y(9);
         TS1 = Y(10); TS2 = Y(11); %hypothetical components involved in the conversion of temperature signal
         
         Temp = square(t*(2*pi)/pt,(pw/pt)*100)*Tempint/2+Tcoreavg;
         
         kin0 = 0.13; kin1 = 0.65; kout1 = 0.65; %from STEP1
         kin2 = 8; Kin2 = 20; n = 3; kout2 = 0.45;
         
         dTS1dt = kin0*Tcoreavg + kin1*(Temp-Tcoreavg) - kout1*TS1;
         dTS2dt = kin2*TS1^n/(Kin2^n+TS1^n)-kout2*TS2;
    
         ks = entrainparam(1); Ks = entrainparam(2);
         khsf = entrainparam(3); Khsf = entrainparam(4);
         
         HSF1tot = hsparam(1); vact = hsparam(2); vina = hsparam(3);
         vbhsp = hsparam(4); vdhsp = hsparam(5);
         
         dactHSF1dt = vact*(HSF1tot-actHSF1)*(1+ks*TS2/(Ks+TS2))/indHSP-vina*actHSF1;
         dindHSPdt = vbhsp*actHSF1/indHSP-vdhsp*indHSP;

         %PCG dynamics and parameter values (Mavroudis et al, 2014)
         v1b = 9; c = 0.01; k1b = 1; k1i = 0.56; p = 8; k1d = 0.12; kc = 0; DRN = 1; 
         k2b = 0.3; q = 2; k2d = 0.05; k2t = 0.24; 
         k3t = 0.02; k3d = 0.12;
         v4b = 3.6; k4b = 2.16; r = 3; k4d = 0.75;
         k5b = 0.24; k5d = 0.06; k5t = 0.45; 
         k6t = 0.06; k6d = 0.12; k6a = 0.09;
         k7a = 0.003; k7d = 0.09;
         
         dPerCry_mRNAdt = v1b*(ClockBmal1+c)*(1+khsf*actHSF1/(Khsf+actHSF1+indHSP/ClockBmal1))/(k1b*(1+(PerCry_nuc/k1i)^p+ClockBmal1+c))-k1d*PerCry_mRNA+kc*DRN/ClockBmal1;
         dPerCrydt = k2b*PerCry_mRNA^q-k2d*PerCry-k2t*PerCry+k3t*PerCry_nuc;
         dPerCry_nucdt = k2t*PerCry-k3t*PerCry_nuc-k3d*PerCry_nuc;
         dBmal1_mRNAdt = v4b*PerCry_nuc^r/(k4b^r+PerCry_nuc^r)-k4d*Bmal1_mRNA;
         dBmal1dt = k5b*Bmal1_mRNA-k5d*Bmal1-k5t*Bmal1+k6t*Bmal1_nuc;
         dBmal1_nucdt = k5t*Bmal1-k6t*Bmal1_nuc-k6d*Bmal1_nuc-k6a*Bmal1_nuc+k7a*ClockBmal1;
         dClockBmal1dt = k6a*Bmal1_nuc-k7a*ClockBmal1-k7d*ClockBmal1;

         dYdt = [dPerCry_mRNAdt;dPerCrydt;dPerCry_nucdt;dBmal1_mRNAdt;dBmal1dt;dBmal1_nucdt;dClockBmal1dt;dactHSF1dt;dindHSPdt;dTS1dt;dTS2dt];
end