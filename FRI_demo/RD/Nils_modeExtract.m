function [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT, Cs_simple] = Nils_modeExtract(x, M, Nr, sigma_s, clwin )
% [ m_SR_Cl,m_SR_MB,m_LCR_Cl, m_LCR_MB, STFT] = Nils_modeExtract(x, M, Nr, sigma_s, clwin )
%
%  Nils mode extraction method
%
% input:
% x : signal
% M : number of frequency bins
% Nr: number of components to extract
% sigma_s : (optional, default 0.09) window width
% clwin : (optional, default =10) frequency clearing window
%
% output (estimated modes):
%
% m_SR_Cl  : Classical method
% m_SR_MB  : 
% m_LCR_Cl : 
% m_LCR_MB : 
% STFT     : STFT matrix
%
%


N = length(x);

if ~exist('sigma_s', 'var')
   sigma_s = 0.09;  
end

if ~exist('clwin', 'var')
   clwin = 10;  
end

if ~exist('create_gaussian_window')
    addpath('./Nils');
end


Nfft = M;


        
[g, Lh] = create_gaussian_window(N, Nfft, sigma_s);
[STFT, omega, ~, QM, ~, tau] = FM_operators(x, N, Nfft, g, Lh, sigma_s);

%% [8] extraction ridge simple (classique ) [R. Carmona,  W. Hwang,  and B. Torresani,  ?Characterization of signals  by the ridges  of their wavelet  transforms,?IEEETransactions on Signal Processing, vol. 45, no. 10, pp. 2586?2590, Oct 1997]
   [Cs_simple] = exridge_mult(STFT, Nr, 0, 0, clwin);

%% reconstruction (m_SR_Cl : simple reconstruction [8], Linear Chirp Reconstruct  (LCR)
   [m_SR_Cl, m_LCR_Cl, IF_Cl] = R1_MR_and_LCR_grid(STFT, QM, Cs_simple, g, Lh, sigma_s, Nr, Nfft, N);

    [Cs_VFB_MB] = VFB_MB_exridge_MCS(STFT, sigma_s, QM, 2, Nr);
    [m_SR_MB, m_LCR_MB, IF_MB] = R1_MR_and_LCR_grid(STFT, QM, Cs_VFB_MB, g, Lh, sigma_s, Nr, Nfft, N);




end

