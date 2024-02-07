clear all
close all
clc

folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'tools']));


MCrep = 100;          % number of MC realization (20-30 are sufficient to reproduce the figure behaviour)

%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 20;        %% analysis window size (in bin)
%% Define signal x0
Ncomp = 3;                          %% number of components
N     = 500;                        %% signal length
X = zeros(N,Ncomp);
% 
X(:,1) = (fmconst(N, 0.1));
X(:,2) = (fmlin(N,0.13,0.3));
X(:,3) = (fmsin(N,0.3,0.45,320,1,0.3,+1));

x0 = sum(X,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 1 SNR 10dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfr = zeros(M,N,MCrep);
Spect = zeros(M/2,N,MCrep);
SNR = 10;
for it = 1:MCrep   %% iterations
    % Add noise
    x = sigmerge(x0, randn(size(x0)), SNR);
    [tfr(:,:,it),stfr] = tfrvsgab2(x, M, L); %% compute SST
    Spect(:,:,it) = abs(stfr(1:M/2,:)).^2;
end

save('sst_SNR10.mat','tfr','Spect')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 1 SNR 0dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfr = zeros(M,N,MCrep);
Spect = zeros(M/2,N,MCrep);
SNR = 0;
for it = 1:MCrep   %% iterations
    % Add noise
    x = sigmerge(x0, randn(size(x0)), SNR);
    [tfr(:,:,it),stfr] = tfrvsgab2(x, M, L); %% compute SST
    Spect(:,:,it) = abs(stfr(1:M/2,:)).^2;
end

save('sst_SNR0.mat','tfr','Spect')



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig2
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfr = zeros(M,N,MCrep);
Spect = zeros(M/2,N,MCrep);
snr_range = [-20 20]; % SNR range to compute
SNRt = snr_range(1):4:snr_range(2);

for indsnr = 1:length(SNRt)
    SNRi = SNRt(indsnr);
    for it = 1:MCrep   %% iterations
        % Add noise
        x = sigmerge(x0, randn(size(x0)), SNRi);
        [tfr(:,:,it,indsnr),stfr] = tfrvsgab2(x, M, L); %% compute SST
        Spect(:,:,it,indsnr) = abs(stfr(1:M/2,:)).^2;
    end
end

save('sst_Fig2tfr.mat','tfr','-v7.3')
save('sst_Fig2Spect.mat','Spect','-v7.3')


