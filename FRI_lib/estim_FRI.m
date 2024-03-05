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
% ifia       : boolean for amplitude estimation
% Oracle     : bolean to use or not the ground truth IF
% tgt        : ground truth IF used for Oracle


% OUTPUT:
% Estim_IF    : Estimated ridges positions
% Estim_A     : Estimated amplitudes
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024



%% Algorithm
[M,N] = size(Spect);                       % Number of time and frequency bins
mv = -M0:M0;                               % we keep 2*M0+1 Fourier coefficients
[~,modF]=max(F);                           % Mode storing
Estim_IF = zeros(N,NbDiracs);              % To store IF estimates
Estim_A = zeros(N,NbDiracs);               % To store amplitude estimates

phihat = (1/(M))*fft(F); % division by LenF (size of F) due to a different definition of the fourier transform
% see the paper 'Super resolved time of flight imaging via FRI sampling theory' for details
phiH(1:2*M0+1) = [phihat(end-M0+1:end);phihat(1:M0+1)];
V = exp((1i * 2 * pi * kron((1:M)',mv)) ./ M);
Dphi = diag(phiH(1:2*M0+1));

%% Estimation
for tt = 1:N
    m = Spect(:,tt); % Current time slice
    y =  Dphi \ (pinv(V) * m);

    %% Annihilating filter method
    h = anF(Method,y,NbDiracs);
    
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

