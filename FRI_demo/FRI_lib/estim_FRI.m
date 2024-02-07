function [Estim_IF,Estim_A] = estim_FRI(Spect,NbDiracs,F,M0,Method,ifia,Oracle,tgt)
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

phihat = (1/(M))*fft(F); % division par LenF (la taille de F) pour obtenir une approximation des coefficients de Fourier
phiH(1:2*M0+1) = [phihat(end-M0+1:end);phihat(1:M0+1)];
V = exp((1i * 2 * pi * kron((1:M)',mv)) ./ M);
Dphi = diag(phiH(1:2*M0+1));

% Estimation
for tt = 1:N
    m = Spect(:,tt); % Load current signal
    y =  Dphi \ (pinv(V) * m);

    %% Annihilating filter method
    h = anF(Method,y,NbDiracs,(length(y)-1)/2,0);
    
    %% Locations retrieval
    [Estim_IF(tt,:)] = FRI_EstimLoc(h,NbDiracs,M);    
    Estim_IF(tt,:) = sort(mod(Estim_IF(tt,:) - (M-modF) -1, M));
    Estim_IF(tt,:) = max(min(Estim_IF(tt,:),M),1);


    %% Estimate amplitude
    if (ifia && ~Oracle)
        [Estim_A(tt,:)] = FRI_EstimWei(mv,Estim_IF(tt,:),NbDiracs,M,y);
    elseif (ifia && Oracle)
        [Estim_A(tt,:)] = FRI_EstimWei(mv,tgt(tt,:),NbDiracs,M,y);
    end

end
