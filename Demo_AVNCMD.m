% ------------------------ Demo for AVNCMD --------------------------------
%
% This is a simple example to test the AVNCMD algorithm 
% -- Adaptive Variational Nonlinear Chirp Mode Decomposition
%
% Author: Hao Liang 
%
% Last modified by: 21/09/19
%

clc; clear; close all

T = 1;         % time duration
fs = 2000;     % sample frequency
t = 0:1/fs:T;  % time variables

% Instantaneous amplitudes (IAs) and instantaneous frequencies (IFs)
a1 = exp(-0.03*t); a2 = exp(-0.06*t);
f_1t = 25 + 8*t - 3*t.^2 + 0.4*t.^3;
f_2t = 40 + 16*t - 6*t.^2 + 0.4*t.^3;

% A two-component simulated nonlinear chirp signal (NCS)
g_1t = a1.*cos(2*pi*(0.8 + 25*t + 4*t.^2 - 1*t.^3 + 0.1*t.^4));
g_2t = a2.*cos(2*pi*(1 + 40*t + 8*t.^2 -2*t.^3+0.1*t.^4));
g = g_1t + g_2t;


%% AVNCMD
beta = 1e-6;  % filter parameter
tol = 1e-5;   % tolerance of convergence criterion
iniIF = [28*ones(1,length(t));48*ones(1,length(t))];  % initial IFs 
         
% Start AVNCMD algorithm
tic;
[estIF, estIA, estMode] = BVNCMD(g, fs, iniIF, beta, tol);
toc;

%% Relative errors of the estimated modes and IFs
IF1_RE =  norm(estIF(1,:,end)-f_1t)/norm(f_1t)
IF2_RE =  norm(estIF(2,:,end)-f_2t)/norm(f_2t)
Mode1_RE =  norm(estMode(1,:,end) - g_1t)/norm(g_1t)
Mode2_RE =  norm(estMode(2,:,end) - g_2t)/norm(g_2t)


%% Show the estimated IF
figure
plot(t,[f_1t;f_2t],'b','linewidth',2); hold on  % true IFs
plot(t,estIF(:,:,end),'r','linewidth',2)  % estimated IFs
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',24,'FontName','Times New Roman');
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'YDir','normal','FontName','Times New Roman')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	
ylim([20 55])
