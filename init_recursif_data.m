function [F,a,b] = init_recursif_data(M,L,k)
% [F,a,b] = init_recursif_data(M,L,k)
%
%  Precompute the squared modulus of 
%  the Fourier transform of the recursive analysis window
% 
%
% input:
% M: number of frequency bins
% L: window time spread 
% k: filter order 
%
% output:
% F: g
% a,b: filter recursive coef

m0 = 0;   %% le bin frequentiel central (ridge)
m_range = (-(M/4)+1:(M/4));

omega0 = 2 * pi * m0/M;
omega = 2*pi* m_range/M;
T = L;
F = (1+((omega-omega0)*T).^2).^(-k);
F = transpose(F/(M/2));


%% module
% omega0 = 2*pi * (M/4)/M;
% omega = -(M/4)+1:(M/4);
% T = L;%L/(M/2);
% F = (1+((omega-omega0)*T).^2).^(-k);
% F = F./max(F);


%% Recursive_STFT
mm = 1:M/2;
lambda = (1+mm-1)/M;
pTs    = -1.0/L + 1i*2*pi.*lambda;
alpha  = exp(pTs);
[a,b]  = Gk22(k, L, alpha);
b = flip(transpose(b),1); a = flip(a(:,2:end),2);



