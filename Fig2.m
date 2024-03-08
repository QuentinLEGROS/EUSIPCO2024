clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Fig.2 of the paper 'Estimation of Instantaneous Frequency and 
% Amplitude of Multi-Component Signals using Sparse Modeling of Signal 
% Innovation', by performing a comparison of the IA estimates RMAE.
%  
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024
%
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'FRI_lib']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'Compute_Amplitude_DF']));

snr_range = [-20 20]; % SNR range to compute
MCrep = 100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length


% Amplitude modulated
amp0(:,1)=(1:-0.5/(N-1):0.5);

X(:,1) = amp0(:,1).*(fmlin(N,0.15,0.3));
tf0(:,1) = ((linspace(0.15,0.3,N)).*M); % IF Ground truth

Ncomp = size(X,2);
x0 = sum(X,2);


%% Method parameters
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 10; % variance of the random walk in the temporal model
n_pad = 50; 
detect=0;
div = 4;
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
PneiMask = 5;
alpha = 0.5;
beta = 0.5;
ifia = 0;
Oracle = 0;

%% Compute ground truth
X = real(transpose(X));
tgt = tf0;

%% Calcul normalisation fft
A = 1/(sqrt(2*pi)*L);
C = -1 / (2*L^2);
K = 2 * L * sqrt(2*log(1/10^(-4)));  %% window length in samples
k = (-round(K/2):round(K/2)).^2;
g = A * exp( C * k);
G = sum(abs(fft(g)))/2;


%% Calcul step_r and step_v
sigm_d = round(sqrt((M/2)/(pi*L)));
step_r = 3*sigm_d;
step_v = 3*step_r;


% Parameter of recursive FRI
k=3;                                        % recursive filter order
[Fr,a,b] = init_recursif_data(M,L,k);

%% analysis window
Mm = M/2; % frequency support considered for the method (select M/2 if you want to take half of the frequency)
meanF = 100;
m = -meanF:(Mm)-meanF-1;                    % frequency support of the convolution kernel
F = transpose(Fh(m, M, L ));                % Convolution kernel
F = F./max(F); 

F_sst=compF_SST(M,200);                    % compute data distribution for SST signal

%% Initialization
methods_name = {'EM',...
                'Oracle EM',...
                'PB',...
                'Oracle PB',...
                'Oracle Local',...
                'FRI LS',...
                'FRI TLS (proposed)',...
                'Oracle FRI TLS',...
                'FRI SST',...
                'FRI recursif'
                };
            
            
methods_to_use = [1 3 5 6 7 8 9 10];   % insert here the indices of the methods to compare (names above)

nb_methods = length(methods_to_use);
SNRt = snr_range(1):2:snr_range(2);

MAE_out = zeros(length(SNRt), nb_methods);

%% Compute RQF
for indsnr = 1:length(SNRt)
  fprintf(1, "+ SNR=%d dB \n", SNRt(indsnr));
  SNRi = SNRt(indsnr);

    for ind_met = 1:length(methods_to_use)
        MAE_tmp = zeros(MCrep,1);
        
        for it = 1:MCrep   %% iterations
            clc;
            disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
            disp(strcat(['SNR : ',num2str(indsnr),' / ',num2str(length(SNRt))]));
            disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))
            
            
            % Add noise
            x = sigmerge(x0, randn(size(x0)), SNRi);
            
            switch(methods_to_use(ind_met))
                case 1  %% Proposed EM 
                        [tfr]  = tfrgab2(real(x), M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                        [~,~,tf,Amp]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'Lap',1e-2,step_r,step_v,ifplot,0,1,G,1,L);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 2  %% Oracle EM
                        [tfr]  = tfrgab2(real(x), M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                        [Amp]=Oracle_EM(Spect',Fc,Ncomp,G,round(tf0),M,L);
                        tf = tf0;
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 3  %% Beta divergence
                        alpha  = 0.4;
                        beta   = 0.2;
                        [tfr]  = tfrgab2(real(x), M, L);
                        [tf,Amp] = PB_Oracle(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,0,PneiMask,tf0);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 4 %% Oracle PB
                        alpha  = 0.4;
                        beta   = 0.2;
                        [tfr]  = tfrgab2(real(x), M, L);
                        [tf,Amp] = PB_Oracle(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,1,PneiMask,round(tf0));
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 5  %% Oracle Local
                        [tfr]  = tfrgab2(x, M, L);
                        tf = tf0;
                        [Amp] = Oracle_Amp_DF(tfr,Ncomp,L,round(tf0)); 
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 6  %% FRI LS
                        Method = 2;
                        M0 = 10;
                        [tfr] = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,Amp] = estim_FRI(Spect,Ncomp,F,M0,Method,1,0,0);
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 7  %% FRI TLS
                        Method = 2;
                        M0 = 10;
                        [tfr] = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tf0);
                        [Amp] = Oracle_Amp_DF(tfr,Ncomp,L,round(tf)); 
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 8  %% Oracle FRI TLS
                        Method = 2;
                        M0 = 10;
                        [tfr] = tfrgab2(x, M, L);
                        Spect = abs(tfr(1:M/2,:)).^2;
                        [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tf0);
                        [Amp] = Oracle_Amp_DF(tfr,Ncomp,L,round(tf)); 
                        I = sort_IF(tfr,tf,Pnei,X,M,L,Ncomp,N,0);
                case 9 %% FRI - SST
                       Method = 2;
                       M0 = 10;
                       [~,stfr, ~,~,~] = tfrvsgab2(real(x), M, L);
                       Spect = abs(stfr(1:M/2,:)).^2;
                       [tf,~] = estim_FRI(Spect,Ncomp,F_sst,M0,Method,ifia,Oracle,tf0);
                       [Amp,x_hat] = Oracle_Amp_DF_stfr(stfr,Ncomp,L,round(tf)); 
                       [I,~] = match_components(X, x_hat); 
                case 10 %% Recursive FRI
                       Method = 2;
                       M0 = 10;
                       tf = RecursiveFRI(x,k,L,n_pad,M,N,a,b,Ncomp,Fr,M0,Method);
                       [tfr] = tfrgab2(x, M, L); %% compute SST
                       [Amp,x_hat] = Oracle_Amp_DF(tfr,Ncomp,L,round(tf)+1);
                       [I,~] = match_components(X, x_hat); 
            end  %% switch

            Amp = Amp(:,I);
            
            % MAE
            MAE_tmp(it) = sum(sum(abs(Amp(n_pad:end-n_pad,:) - amp0(n_pad:end-n_pad,:))))/(N*Ncomp);
        end    %% for methods
        MAE_out(indsnr, ind_met) = mean(MAE_tmp(~isnan(MAE_tmp)));
    end  %% methods
 
end %% snrs

% Normalization
MAE_out = MAE_out ./(N*Ncomp);


%% Plot
cols         = {'k-x' 'b-x' 'g-x' 'r-x' 'k-o' 'b-s' 'g-.' 'r-.' 'b-v'  'b--' 'b*'};
leg = {};

figure(1)
iii = 0;
for ind_met =  1:nb_methods
 if ind_met == 1
  hold off
 else
  hold on
 end

 h(ind_met) = plot(SNRt, squeeze(MAE_out(:,ind_met)), cols{ind_met});
 xlabel('SNR (dB)', 'FontSize', 15)
 ylabel('RMAE', 'FontSize', 15)
 leg{ind_met} = methods_name{methods_to_use(ind_met)};
end
legend(h, leg, 'location', 'NorthWest', 'FontSize', 8)
grid
axis square
set(gca, 'YScale', 'log')
