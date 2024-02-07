function [Estim_IF,Estim_A] = estim_FRI_modulation(Spect,NbDiracs,F,M0,Method,L,tgt)
% [Estim_Spect,Estim_IF,Estim_A] = estim_FRI(Spect,NbDiracs,F,M0)
%
% Main algorithm - estimate the modes IFs and amplitudes using the FRI
% principle. The presence of noise is handled using state-of-the-art
% algorithms
%
% INPUT:
% Spect      : Cut spectrogram - only the first half frequencies 
% NbDiracs   : Number of components
% F          : Filtering function
% M0         : Number of kept frequency components
% choice     : For additional denoising step (1:without | 2:with)
% Method     : For replacing the Prony method in the presence of noise (1:without | 2:TLSD | 3:Cadzow denoising)
%
% OUTPUT:
% Estim_IF    : Estimated ridges positions
% Estim_A     : Estimated amplitudes
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 12-may-2021



%% Algorithm
[M,N] = size(Spect);                       % Number of time and frequency bins
mv = -M0:M0;                               % we keep 2*M0+1 Fourier coefficients
[~,modF]=max(F);                           % Mode storing
Estim_IF = zeros(N,NbDiracs);              % To store IF estimates
Estim_A = zeros(N,NbDiracs);               % To store amplitude estimates

mm = (1:M);
phihat = 1/M*fft(F); % division par LenF (la taille de F) pour obtenir une approximation des coefficients de Fourier
phiH(1:2*M0+1) = [phihat(end-M0+1:end);phihat(1:M0+1)];
V = exp((1i * 2 * pi * kron((1:M)',mv)) ./ M);
Dphi = diag(phiH(1:2*M0+1));


% First time bin
tbin = 10;
m = Spect(:,tbin); % Load current signal
y =  Dphi \ (pinv(V) * m);
h = anF(Method,y,NbDiracs,(length(y)-1)/2,0);
[Estim_IF(tbin,:)] = FRI_EstimLoc(h,NbDiracs,M,modF); 
Estim_IF(tbin,:) = sort(mod(Estim_IF(tbin,:) - (M-modF) -1, M))-1;
tf0 = round(Estim_IF(tbin,:));
Estim_A(tbin,:) = FRI_EstimWei(mv,Estim_IF(tbin,:),NbDiracs,M,y,modF);

% Next time bins
for tt = 2:N
    m = Spect(:,tt); % Load current signal
    y =  Dphi \ (pinv(V) * m);
    
    %% Annihilating filter method
    h = anF(Method,y,NbDiracs,(length(y)-1)/2,0);
    
    %% Locations retrieval
    [Estim_IF(tt,:)] = FRI_EstimLoc(h,NbDiracs,M,modF);    
    Estim_IF(tt,:) = sort(mod(Estim_IF(tt,:) - (M-modF) -1, M))-1;
    % tf = round(Estim_IF(tt,:));

    %% Frequency modulation
    tfdiff = (tf-tf0)./M;
    % [tf0 tf]
    for Nc = 1:NbDiracs
        [Fchirp(:,Nc),A] = (stft_chirp(mm, 1, M, L, tf0/M, tfdiff(Nc)));
        Fchirp= Fchirp.^2;   %% square modulus
    end
    tf0 = tf;

    phihat = 1/M*fft(F); % division par LenF (la taille de F) pour obtenir une approximation des coefficients de Fourier
    phiH(1:2*M0+1) = [phihat(end-M0+1:end);phihat(1:M0+1)];

    Dphi = diag(phiH(1:2*M0+1));
    y =  Dphi \ (pinv(V) * m);

% figure(10)
% hold on
% plot(m)
% plot(F, 'g-.')
% plot(Fchirp, 'r--')
% hold off
% 
% 
% A


    %% Estimate amplitude
    [Estim_A(tt,:)] = FRI_EstimWei(mv,Estim_IF(tt,:),NbDiracs,M,y,modF);
        % [Estim_A(tt,:)] = FRI_EstimWei(mv,tgt(tt,:),NbDiracs,M,y,modF);

% close all
end




end
