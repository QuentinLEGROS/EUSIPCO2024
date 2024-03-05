% clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Table.1 of the paper 'Estimation of Instantaneous Frequency and 
% Amplitude of Multi-Component Signals using Sparse Modeling of Signal 
% Innovation', by comparing the execution time of the compared methods.
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'Brevdo']));
addpath(strcat([folder 'mfiles']));
addpath(strcat([folder 'FRI_lib']));
addpath(strcat([folder 'RecursiveSTFT']));
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'PseudoBay']));
addpath(strcat([folder 'Modul_EM']));
addpath(strcat([folder 'RD']));

snr_range = [-20 20]; % SNR range to compute
MCrep =10;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
% M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
N     = 500;                        %% signal length

clear X tf0

% First component - sinusoidal
X(:,1) = (fmconst(N, 0.1));

% Second component - linear chirp
X(:,2) = (fmlin(N,0.13,0.3));

% Third component - FM chirp
[X(:,3),tf0(:,3)] = (fmsin(N,0.3,0.45,320,1,0.3,+1));

x0 = sum(X,2);
Ncomp = size(X,2);                          %% number of components
X = real(transpose(X));

%% Method parameters
Pnei = 10; % the output mask is of size 2*Pnei+1 -> extract only the ridge position
ifplot =  0; % plot intermediary estimation of each ridge using the proposed approach
ds    = 5; % variance of the random walk in the temporal model
alpha = 0.5;
beta = 0.5;
detect = 0;
ifia = 0;
Oracle = 0;
tgt = tf0;
n_pad = 50;

%% Calcul normalisation fft
A = 1/(sqrt(2*pi)*L);
C = -1 / (2*L^2);
K = 2 * L * sqrt(2*log(1/10^(-4)));  %% window length in samples
k = (-round(K/2):round(K/2)).^2;
g = A * exp( C * k);
G = sum(abs(fft(g)))/2;







%% Initialization
methods_name = {'Brevdo',...
                'PB',...
                'RD',...
                'EM', ... 
                'FRI TLS',...
                'Recursive FRI'
                };
            
methods_to_use = 6%[1 2 3 5 6];   % insert here the indices of the methods to compare (names above)
nb_methods = length(methods_to_use);
mm = [500 500 1000 2000];
% ctime = zeros(MCrep,1);
Ctime = zeros(length(methods_to_use),length(mm));


%% Compute RMSE

for ind_met = 1:length(methods_to_use)

    for ind_m = 1:length(mm)
        M = mm(ind_m);

    for it = 1:MCrep   %% iterations
        clc;
        disp(strcat(['Method :', methods_name{methods_to_use(ind_met)}]));
        disp(strcat(['M : ', num2str(mm(ind_m))]))
        disp(strcat(['Iter : ', num2str(it),' / ',num2str(MCrep)]))


        % Add noise
        x = sigmerge(x0, randn(size(x0)), 10);
        [tfr]  = tfrgab2(x, M, L);
        Spect = abs(tfr(1:M/2,:)).^2;
        
        switch(methods_to_use(ind_met))
            case 1  %%Brevdo STFT
                tic
                [~,mask,tf] = Brevdo_modeExtract(tfr, L, Ncomp, Pnei);
                ctime(it) = toc;
            

             case 2  %% PB
                beta  = 0.5; % beta divergence hyperparameter  ||  POSITIVE AND DIFFERENT TO 0
                div   = 2;   % 3 = Renyi
                PneiMask = 5;
                tic
                [mask,tf] = PB_Oracle(tfr,Ncomp, M, L, div, beta, alpha, ds, Pnei, ifplot,detect,0,PneiMask);
                ctime(it) = toc;

            case 3  %% RD
                 Nr = Ncomp;
                sigma_s = 0.09;
                clwin = 10;
                tic
                [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT, Cs_simple] = Nils_modeExtract(x, M, Nr, sigma_s, clwin );
                ctime(it) = toc;
              
             case 4 %% EM TV
                 %% Calcul step_r and step_v
                sigm_d = round(sqrt((M/2)/(pi*L)));
                step_r = 3*sigm_d;
                step_v = 3*step_r;
                [Fc]=comp_Fc(M,L);Fc = Fc + eps;     %% Data distribution
                tic
                [~,~,tf,~]=Mod_Estim_W_EM_multi(Spect',Fc,Ncomp,1,'TV',1e-4,step_r,step_v,ifplot,0,1,G,1,L);
                ctime(it) = toc;

             case 5 % FRI TLS
                Mm = M/2; % frequency support considered for the method (select M/2 if you want to take half of the frequency)
                meanF = 125;
                m = -meanF:(Mm)-meanF-1;                    % frequency support of the convolution kernel
                F = transpose(Fh(m, M, L ));                % Convolution kernel
                Method = 2;
                M0 = 20;
                [tfr] = tfrgab2(x, M, L); %% compute SST
                Spect = abs(tfr(1:M/2,:)).^2;
                tic
                [tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,ifia,Oracle,tgt);
                ctime(it) = toc;
            case 6 % Recursive FRI
                k=3;                                        % recursive filter order
                [Fr,a,b] = init_recursif_data(M,L,k);
                Method = 2;
                M0 = 10;
                tic
                tf = RecursiveFRI(x,k,L,n_pad,M,N,a,b,Ncomp,Fr,M0,Method);
                % n0 = ((k-1)*L);
                % xt = [x;zeros(n_pad+1,1)]; % padding
                % [tfr, nfreqs] = recursive_stft(xt, k, L, 1, M, M); % calcul stft récursive
                % Spect = abs(tfr(1:M/2,:)).^2;
                % tic
                % parfor n = k:length(xt)
                %      [tf(n,:)] = estim_FRI_recursif(Spect(:,n),Ncomp,Fr,M0,Method);
                % end
                ctime(it) = toc;

        end  %% switch
    end %% repetitions


    Ctime(ind_met,ind_m) = mean(ctime);
    end
end  %% methods


for ind_met = 1:length(methods_to_use)
    for ind_m = 1:length(mm)
        disp(strcat([methods_name(ind_met),', M = ',num2str(mm(ind_m)),' : ',num2str(Ctime(ind_met,ind_m)) ]))
    end
end
 



