clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Fig.2 of the paper ESTIMATION OF OVERLAPPING MODES USING SPARSE MODELING OF SIGNALINNOVATION
%  Comparison of the RMSE in the non-overlapping case
%
%  Authors : Q.Legros (quentin.legros@telecom-paris.fr) and D.Fourer
%  Date    : 16-feb-2021
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'FRI_lib']));
addpath(strcat([folder 'PseudoBay']));


snr_range = [-20 20]; % SNR range to compute

%% Time-frequency representation parameters
N     = 500;       %% signal length
M     = 500;       %% Number of frequential bin
L     = 20;        %% analysis window size (in bin)
npad  = 50;
%% Define signal x0
Amp(:,1) = 2*ones(N,1);%linspace(1,5,N);
Amp(:,2) = linspace(2,1,N);
% Amp(:,3) = 1*ones(N,1);%linspace(2,1,N);

X(:,1) = Amp(:,1).*(fmconst(N, 0.1));
tf0(:,1) = 0.1*ones(N,1);
% X(:,2) = Amp(:,2).*(fmlin(N,0.15,0.4));
% tf0(:,2) = linspace(0.25,0.4,N);
X(:,2) = Amp(:,2).*(fmconst(N, 0.4));
tf0(:,2) = 0.4*ones(N,1);
% X(:,3) = (fmsin(N,0.3,0.45,320,1,0.3,+1));

Ncomp = size(X,2);                          %% number of components


x0 = sum(X,2);
X = transpose(X);

%% Method parameters
M0 = 10;                                    % Frequency truncation - to avoid infinite sum
[F_mat]=comp_Fc(M,L);
F = F_mat(:,200); % F = F./max(F);                         % truncate to data size


%% Compute ground truth
[tfr]  = tfrgab2(x0, M, L);
spect=(abs(tfr(1:round(M/2),:))).^2;

tf0 = sort(tf0,2);

%% Main

SNR = 10;
% Add noise
x = sigmerge(x0, randn(size(x0)), SNR);

Method = 2;
[tfr] = tfrgab2(x, M, L); %% compute SST
Spect = abs(tfr(1:M/2,:)).^2;
[tf,ia] = estim_FRI(Spect,Ncomp,F,M0,Method);

figure(1)
subplot(2,1,1)
imagesc(flipud(spect))
subplot(2,1,2)
plot(tf)
ylim([0 250])

figure(2)
title('amplitudes')
subplot(2,1,1)
hold on
plot(ia(npad:end-npad,1))
plot(Amp(npad:end-npad,1))
hold off
legend('estimation','Ground truth')
subplot(2,1,2)
hold on
plot(ia(npad:end-npad,2))
plot(Amp(npad:end-npad,2))
hold off
legend('estimation','Ground truth')

