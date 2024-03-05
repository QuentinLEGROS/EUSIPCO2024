function [Estim_IF] = estim_FRI_recursif(Spect,NbDiracs,F,M0,Method)
% [Estim_Spect,Estim_IF,Estim_A] = estim_FRI(Spect,NbDiracs,F,M0)
%
% Main algorithm - estimate the modes IFs and amplitudes using the FRI
% principle. The presence of noise is handled using state-of-the-art
% algorithms.
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
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024
%

%% precomputation
[M,~] = size(Spect);                       % Number of time and frequency bins
mv = -M0:M0;                               % we keep 2*M0+1 Fourier coefficients
[~,modF]=max(F);                           % Mode storing

phihat = 1/M*fft(F); % division par LenF (la taille de F) pour obtenir une approximation des coefficients de Fourier
phiH(1:2*M0+1) = [phihat(end-M0+1:end);phihat(1:M0+1)];
V = exp((1i * 2 * pi * kron((1:M)',mv)) ./ M);
Dphi = diag(phiH(1:2*M0+1));



%% Main algorithm
y =  Dphi \ (pinv(V) * Spect);

%% Annihilating filter method
h = anF(Method,y,NbDiracs);

%% Locations retrieval
[Estim_IF] = FRI_EstimLoc(h,NbDiracs,M);    
Estim_IF = sort(mod(Estim_IF - (M-modF) -1, M));
Estim_IF = max(min(Estim_IF,M),1);

end


