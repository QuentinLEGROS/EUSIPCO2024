function tf = RecursiveFRI(x,k,L,n_pad,M,N,a,b,Ncomp,Fr,M0,Method)
% tf = RecursiveFRI(x,k,L,n_pad,M,N,a,b,Ncomp,Fr,M0,Method)
%
% Main algorithm - estimate the modes IFs and amplitudes using the Recursive FRI
% method. The STFT is computed through a recursive implementation before
% applying the FRI TLS algorithm.
%
% INPUT:
% x           : signal
% k           : degree for recursive implementation of the STFT
% L           : time spread parameter of the analysis window
% n_pad       :   padding for the STFT computation
% M           : number of frequency bin
% N           : signal length
% a           : parameter for the recursive implementation of the STFT
% b           : parameter for the recursive implementation of the STFT
% Ncomp       : Number of component
% Fr          : analysis window
% M0          : Number of kept frequency components
% Method      : For replacing the Prony method in the presence of noise (1:without | 2:TLSD | 3:Cadzow denoising)
%
% OUTPUT:
% tf          : Estimated ridges positions
%
% Author      : Q.Legros (quentin.legros@univ-orleans.fr) and D. Fourer
% Date        : 1-mar-2024


n0 = ((k-1)*L);
xt = [zeros(k-1,1);x;zeros(n_pad+1,1)]; % padding
tfrp    = zeros(M/2, N+n_pad+1);
Spect = zeros(M/2, N+k);
tf = zeros(N+n_pad+k,Ncomp);


% Calcul IF via FRI
for n = n_pad-n0:length(xt)
    tfrp(:,n+1) = transpose(sum(b.*xt(n-k+1:n),1)) - sum(a.*tfrp(:,n-k+1:n),2);
    Spect(:,n-k+1) = abs(tfrp(:,n+1)).^2;
    [tf(n,:)] = estim_FRI_recursif(Spect(:,n-k+1),Ncomp,Fr,M0,Method);
end

% tf = [tf(end-N-k-(k-1):end-k-(k-1)-1,:)];
tf = [tf(end-N-1:end-k+1,:)];
tf = max(min(tf,M/2),1);