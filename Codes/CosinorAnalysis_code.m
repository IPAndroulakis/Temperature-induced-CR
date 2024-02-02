%% Cosinor Analysis
clear all; clc; close all;

t_exp = [12 21 24 28 44 48 53 68 72 77 92 96 101 108 116 120 126 132 140 144];
PER2_exp = [0.98 1.19 1.01 0.8 1.17 1.05 0.8 1.16 1.04 0.79 1.17 1.01 0.79 1.03 1.18 1.01 0.79 1.01 1.19 1.01]; %Fig.1C (Saini et al, 2012), "data for cosinor analysis.doc"

x = t_exp;
y = PER2_exp;
yub = max(y);
ylb = min(y);
yrange = (yub-ylb);                                 % Range of y

phase = x(2)-6;                                     % Estimate phase
ym = mean(y);                                       % Estimate offset
%{
fit = @(b,x)  b(1).*cos(2*pi*(x-b(2))./b(3)) + b(4);           % Function to fit: b = [amp,phase,period,mean]
fcn = @(b) sum((fit(b,x) - y).^2);                             % Least-Squares cost function
s = fminsearch(fcn, [yrange; phase; period; ym])               % Minimise Least-Squares
%}
fit = @(b,x)  b(1).*cos(2*pi*(x-b(2))./24) + b(3);           % Function to fit: b = [amp,phase,mean]
fcn = @(b) sum((fit(b,x) - y).^2);                           % Least-Squares cost function
[s,fval,exitflag] = fminsearch(fcn, [yrange; phase; ym])     % Minimise Least-Squares
%xfit = linspace(min(x),max(x),500);
xfit = linspace(0,max(x),500);

figure(1)
plot(x,y,'ko', xfit,fit(s,xfit),'r-');
hold on

%% Grab points from Cosine curve
t_exp_cosinor = 0:0.5:24;
PER2_exp_cosinor = s(1).*cos(2*pi*(t_exp_cosinor-s(2))./24) + s(3);  
figure(1)
plot(t_exp_cosinor,PER2_exp_cosinor,'b.');
hold on
xlabel('Time (h)'); xticks(0:12:144)
ylabel('PER2 concentration (a.u.)')
legend('Exp data','Fitted cosine curve','Grabbed data points');
set(gca,'FontSize',16);

data_exp_cosinor = [t_exp_cosinor', PER2_exp_cosinor'];
writematrix(data_exp_cosinor,'data_exp_cosinor.xlsx','WriteMode','overwrite');
