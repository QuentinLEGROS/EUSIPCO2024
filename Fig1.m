clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Fig.1 of the paper 'Estimation of Instantaneous Frequency and 
% Amplitude of Multi-Component Signals using Sparse Modeling of Signal 
% Innovation', by performing a comparison of the IF estimates RMSE.
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'RecursiveSTFT']));
addpath(strcat([folder 'FRI_lib']));
addpath(strcat([folder 'RD']));

snr_range = [-20 20]; % SNR range to compute
MCrep = 100; %100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

% First component - sinusoidal
X(:,1) = (fmconst(N, 0.1));
tf0(:,1)= (0.1.*ones(N,1) .* M)+1; % IF Ground truth

% Second component - linear chirp
X(:,2) = (fmlin(N,0.13,0.3));
tf0(:,2) = ((linspace(0.13,0.3,N)).*M) +1; % IF Ground truth

% Third component - FM chirp
[X(:,3),tf0(:,3)] = (fmsin(N,0.3,0.45,320,1,0.3,+1));
tf0(:,3) = (tf0(:,3) .* M) + 1; % IF Ground truth

% Number of component
Ncomp = size(X,2);
x0 = sum(X,2);
X = transpose(X);


%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifia = 0;
Oracle = 0;
tgt = tf0;

% window analysis for FRI methods
Mm = M/2; % frequency support considered for the method (select M/2 if you want to take half of the frequency)
meanF = 125;
m = -meanF:(Mm)-meanF-1;                    % frequency support of the convolution kernel
F = transpose(Fh(m, M, L ));                % Convolution kernel

% window analysis for FRI SST
F_sst=compF_SST(M,200);                    % compute data distribution for SST signal
n_pad = 50;                                % padding for comparison


% Parameter of PB methods
alpha  = 0.5;
beta   = 0.5;
ds     = 2; % variance of the random walk in the temporal model
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach

% Parameter of recursive FRI
k=3; %3;                                        % recursive filter order
[Fr,a,b] = init_recursif_data(M,L,k);


%% Initialization
methods_name = {'Brevdo [10]',...
                'PB-ABD-\alpha=0.4,\beta=0.2 [21]',...
                'PB-\beta=0.7,Beta-div. [9]',...
                'PB-\alpha=0.3,Renyi-div. [9]',...
                'RD [8]', ...  
                'FRI',...
                'FRI TLS, M0 = 5 (proposed)',...
                'FRI TLS, M0 = 10 (proposed)',...
                'FRI TLS, M0 = 20 (proposed)',...
                'FRI SST (proposed)',...
                'Recursive FRI (proposed)'
                };

methods_to_use = [1 2 3 4 5 6 7 8 9 10 11];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):2:snr_range(2);

%% Compute MAE
L2ErPos_out = zeros(length(SNRt), nb_methods);

for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        L2ErPos_tmp = zeros(length(MCrep));

        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            

            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            
            switch(methods_to_use(ind_met))
               case 1 % Brevdo
                        [tfr] = tfrgab2(x, M, L);       %% compute STFT
                        [~,~,tf] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);
                        tf = tf';   
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
               case 2  %% alpha Beta divergence
                        alpha  = 0.4;
                        beta   = 0.2;
                        div   = 4; 
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
               case 3  %% Beta divergence
                        beta  = 0.7; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                        div   = 2;   % 2 = beta
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 4  %% Reyni divergence
                        alpha = 0.3; % Renyi divergence hyperparameter ||  POSITIVE AND DIFFERENT TO 1
                        div   = 3;
                        [tfr]  = tfrgab2(x, M, L);
                        [~,tf] = pseudoBay(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei,ifplot);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
               case 5  %% simple method
                        Nr = Ncomp;
                        sigma_s = 0.09;
                        clwin = 10;
                        [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT, Cs_simple] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                        tf = Cs_simple';
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 6  %% FRI
                        Method = 1;
                        M0 = 10;
                        [tfr] = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tgt);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 7  %% FRI TLS
                        Method = 2;
                        M0 = 5;
                        [tfr] = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tgt);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 8  %% FRI TLS
                        Method = 2;
                        M0 = 10;
                        [tfr] = tfrgab2(x, M, L); %% compute SST
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tgt);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 9  %% FRI TLS
                        Method = 2;
                        M0 = 20;
                        [tfr] = tfrgab2(x, M, L); %% compute SST
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tgt);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 10 %% FRI - SST
                       Method = 2;
                       M0 = 10;
                       [~,stfr, ~,~,~] = tfrvsgab2(real(x), M, L);
                       Spect = abs(stfr(1:M/2,:)).^2;
                       [tf,~] = estim_FRI(Spect,Ncomp,F_sst,M0,Method,ifia,Oracle,tgt);
                       [tfr] = tfrgab2(x, M, L);
                       I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 11 %% Recursive FRI
                       Method = 2;
                       M0 = 10;
                       tf = RecursiveFRI(x,k,L,n_pad,M,N,a,b,Ncomp,Fr,M0,Method);
                       [tfr] = tfrgab2(x, M, L); %% compute SST
                       I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
            end  %% switch

            % sorting modes for evaluation
            tf = tf(:,I);

            L2ErPos_tmp(it) = sum(sum((tf(n_pad:(end-n_pad),:)-tf0(n_pad:(end-n_pad),:)).^2)./(M*M));
        end    %% repetitions
        L2ErPos_out(indsnr, ind_met) = mean(L2ErPos_tmp(~isnan(L2ErPos_tmp)));
    end  %% methods
 
end %% snrs


%% Plot -- 1
cols         = {'k-x' 'r-o' 'r-v' 'r-x' 'k-^' 'g-o' 'b-v' 'b-x' 'g-x' 'g-^' 'g-+'};
leg = {};

figure(1)
for ind_met =  1:nb_methods
    
 if ind_met == 1
  hold off
 else
  hold on
 end
 
 h(ind_met) = plot(SNRt, 10*log10(squeeze(L2ErPos_out(:,ind_met))), cols{ind_met});
 xlabel('SNR [dB]', 'FontSize', 14)
 ylabel('RMSE [dB]', 'FontSize', 14)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end
legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)
set(gca, 'YScale', 'log')
grid




save('FigIF')



