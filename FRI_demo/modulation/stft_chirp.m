function [ res,A ] = stft_chirp(m, n, M, L, f0, q)
% [ res ] = stft_chirp(m, n, M, L, f0, q)
% 
%  Compute the modulus of the STFT of a linear chirp
%
%
%  INPUT: 
%  m : frequency bin(s) to compute
%  n : time instants to considers
%  M : number of frequency bins to process, fft size (default: length(x))
%  L : time spread of the analysis window
%  f0: initial frequency of the chirp at n=0
%  q : chirp rate
%
%  OUTPUT
%  res : resulting modulus of the STFT
%
%  Author: Dominique Fourer (dominique@fourer.fr)
%  Date : 10-june-2021

N =  length(n);
res = zeros(length(m), N);


lambda2 = (((2*pi*q)^2*L^4)+1);  %% squared modulus of lambda

phi_p = 2*pi*(f0+q*n);           %% instantaneous frequency
omeg = 2*pi*m/M;                 %% frequency bin

A = (lambda2)^(-1/4);
for n0 = 1:N
  res(:,n0) = A .* exp(-L^2*(phi_p(n0)-omeg).^2/(2*lambda2));
end

end

