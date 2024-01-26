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
addpath(strcat([folder 'RecursiveSTFT']));
addpath(strcat([folder 'modulation']));

%% Time-frequency representation parameters
N     = 500;       %% signal length
M     = 500;       %% Number of frequential bin
L     = 20;        %% analysis window size (in bin)
npad  = 50;




%% Define signal x0
Amp(:,1) = 1*ones(N,1);%linspace(1,5,N);
Amp(:,1) = linspace(2,1,N);
% Amp(:,3) = 1*ones(N,1);%linspace(2,1,N);

X(:,1) = Amp(:,1).*(fmconst(N, 0.1));
tf0(:,1) = 0.1*ones(N,1);
% X(:,1) = Amp(:,1).*(fmlin(N,0.15,0.4));
% tf0(:,1) = linspace(0.15,0.4,N);
% X(:,2) = Amp(:,2).*(fmconst(N, 0.4));
% tf0(:,2) = 0.4*ones(N,1);
% X(:,3) = (fmsin(N,0.3,0.45,320,1,0.3,+1));

Ncomp = size(X,2);                          %% number of components
x0 = sum(X,2);
X = transpose(X);



%% Method parameters
M0 = 10;                                    % Frequency truncation - to avoid infinite sum
[F_mat]=comp_Fc(M,L);
F = F_mat(:,200); % F = F./max(F);                         % truncate to data size


% m = -(M/4):(M/4)-1;                       % frequency support of the convolution kernel
% F2 = transpose(Fh(m, M, L ));              % Convolution kernel
% F2 = F2./sum(F2);
% 
% 
% figure
% hold on
% plot(F)
% plot(F2)
% hold off

Method = 1;

%% Compute ground truth
[tfr]  = tfrgab2(x0, M, L);
Spect0=(abs(tfr(1:round(M/2),:))).^2;



%% Recursive_STFT
k=3;                   %% recursive filter order
mi=1; mf=250;          %% frequency range

Nfbins = mf-mi+1;
nfreqs = (mi:mf)/M;
tfr    = zeros(Nfbins, N);



%% Noise
SNR = inf;
x = sigmerge(x0, randn(size(x0)), SNR);

% 
% for m = 1:Nfbins,
%     lambda = (mi+m-1)/M;
%     pTs    = -1.0/L + 1i*2*pi*lambda;
%     alpha  = exp(pTs);
%     [a,b]  = Gk2(k, L, alpha);
%     tfr(m,:) = filter(b,a,x);
% 
%     % imagesc(abs(tfr(1:m,:)).^2)
%     % ylim([1:M/2])
%     % pause(0.05)
% end





%% Main
[tfr] = tfrgab2(x, M, L); %% compute SST
Spect = abs(tfr(1:M/2,:)).^2;
[tf,ia] = estim_FRI(Spect,Ncomp,F,M0,Method,L,tf0(:,1).*M);



%% Plots
figure(1)
subplot(2,1,1)
imagesc((Spect0))
subplot(2,1,2)
plot(tf)
ylim([0 250])

figure(2)
title('amplitudes')
% subplot(2,1,1)
hold on
plot(ia(npad:end-npad,1))
plot(Amp(npad:end-npad,1))
hold off
legend('estimation','Ground truth')
% subplot(2,1,2)
% hold on
% plot(ia(npad:end-npad,2))
% plot(Amp(npad:end-npad,2))
% hold off
% legend('estimation','Ground truth')

