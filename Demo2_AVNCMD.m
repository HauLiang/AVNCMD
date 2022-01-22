% ------------------------ Demo for AVNCMD --------------------------------
%
% This is a new example to test the AVNCMD algorithm 
% -- Adaptive Variational Nonlinear Chirp Mode Decomposition
% -- Advised by a reviewer
%
% Author: Hao Liang 
%
% Last modified by: 22/01/21
%

clc; clear; close all

T = 1;         % time duration
fs = 2000;     % sample frequency
t = 0:1/fs:T;  % time variables

% Instantaneous amplitudes (IAs) and instantaneous frequencies (IFs)
a1 = 1+0.5*cos(2*pi*t); a2 = 1-0.5*cos(2*pi*t);
f_1t = 100 + 300*t + 9*pi*cos(9*pi*t);
f_2t = 400 - 300*t + 9*pi*cos(9*pi*t);

% A two-component simulated nonlinear chirp signal (NCS)
g_1t = a1.*cos(2*pi*(100*t+150*t.^2+sin(9*pi*t)));
g_2t = a2.*cos(2*pi*(400*t-150*t.^2+sin(9*pi*t)));
g = g_1t + g_2t;


%% AVNCMD
beta = 1e-8;  % filter parameter
tol = 1e-5;   % tolerance of convergence criterion
iniIF = [100+300*t;400-300*t];  % initial IFs 

iniIF1 = [100+300*t;400-300*t];
% Start AVNCMD algorithm
tic;
[estIF, estIA, estMode] = AVNCMD(g, fs, iniIF, beta, tol);
toc;


%% Relative errors of the estimated modes and IFs
IF1_RE =  norm(estIF(1,:,end)-f_1t)/norm(f_1t)
IF2_RE =  norm(estIF(2,:,end)-f_2t)/norm(f_2t)
Mode1_RE =  norm(estMode(1,:,end) - g_1t)/norm(g_1t)
Mode2_RE =  norm(estMode(2,:,end) - g_2t)/norm(g_2t)


%% Show the estimated IF
figure; x1 = 0.4; y1 = 200; x2 = 0.6; y2 = 300;
plot(t,[f_1t;f_2t],'b','linewidth',1.5);
hold on;
plot(t,estIF(:,:,end),'r','linewidth',1.5);
set(gcf,'Position',[20 100 640 500]);	
set(gcf,'Color','w'); 
xlabel('Time (s)','FontName','Times New Roman');
ylabel('Frequency (Hz)','FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
ylim([0 500])
set(gca,'FontSize',24)
set(gca,'linewidth',2);
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','k','Linewidth',1);
h1 = axes('position',[0.62 0.22 0.25 0.25]);
axis(h1);
plot(t,[f_1t;f_2t],'b','linewidth',1.5);
hold on;
plot(t,estIF(:,:,end),'r','linewidth',1.5);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12)