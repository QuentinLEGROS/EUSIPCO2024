function [Estim_IF,Estim_A] = estim_FRI(Spect,NbDiracs,F,M0,Method)
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

for tt = 1:N
% for tt = 220:270
    m = Spect(:,tt); % Load current signal
    
    phihat = 1/M*fft(F); % division par LenF (la taille de F) pour obtenir une approximation des coefficients de Fourier
    phiH(1:2*M0+1) = transpose(conj([phihat(end-M0+1:end)' phihat(1:M0+1)']));

    %% Calcul y_m
    % Approximation aux moindres carrés de y selon l'article 'Super resolved time of flight imaging
    % via FRI sampling theory'
    V = exp((1i * 2 * pi * kron((1:M)',mv)) ./ M);
    Dphi = diag(phiH(1:2*M0+1));
%     Dphi = diag(1./(phiH(1:2*M0+1)));
    % y 
    y =  Dphi \ (pinv(V) * m);
    
%     figure(2)
%     plot(abs(pinv(V) * m))
%     pause(0.01)
%     title(strcat(['n = ',num2str(tt)]))
%     
    %% Annihilating filter method

    if Method == 1 % Yulle-Walker system
        h = YW(y);

    elseif Method == 2 % TLSM
        h = AF_TLSM(y,NbDiracs,(length(y)-1)/2,0)';
        
    else
        error('Select Method in [1,3]')
    end

    % Locations retrieval
    [Estim_IF(tt,:)] = FRI_EstimLoc(h,NbDiracs,M,modF);    
    
    % Weight retrvial
    %[Estim_A(tt,:)] = FRI_EstimWei(mv,Estim_IF(tt,:),NbDiracs,M,y);
    
    % Correction according to the data distribution shift (not centered)
    Estim_IF(tt,:) = sort(mod(Estim_IF(tt,:) - (M-modF) -1, M));
    Estim_IF = round(Estim_IF);

    Estim_IF = max(min(Estim_IF,M),1);
    
end


















