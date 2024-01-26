%% Demo (modified version by DF)
clear
close all
clc




folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'modulation']));



%% Time-frequency representation parameters
M       = 500;       %% Number of frequential bin
L       = 50;         %% analysis window size (in bin)

%% Define signal x0
N     = 500;                         %% signal length


f1 = 0.13;                   %% start freq at bin 1
f2 = 0.3;                    %% final freq at bin N
q1 = (f2 - f1)/(N-1)         %% chirp rate (estimated from f2 and f1

%% time instants
n = (1:N)-1;
%% frequencies only from 0 to M/2-1)
m = (1:round(M/2));    %% le bin 1 de matlab equivaut a la frequence 0


% fmlin(N,f1, f2) %% fmlin should be equivalent to the equation below:
x0    = exp(2*pi*1j*(f1 * n + q1 * n.^2/2)); %real
% tf = linspace(f1,f2,N);
% tfdiff = tf(2:end)-tf(1:end-1);
% tfdiff(100)

%% Computation of the STFT
[tfr]  = tfrgab2(x0, M, L);          
spec = abs(tfr(m,:)).^2;   %% spectrogram (squared modulus of the STFT

figure
imagesc(abs(tfr).^2);

%% select a time instant at the middle of the signal
n0 = round(N/2);


figure(1)
plot(m-1, spec(:,n0),'k-')   %% spectrogram of the observed wave sine 

%% locate the ridge at m0 (local maximum)
[~,m0] = max(spec(:,n0));

hold on
%% display the theoretical value of the window Fh (model 1: sinusoid)
plot(m-1, Fh(m-m0,M,L), 'g-.')  

%% Compute the modulus of the STFT of a linear chirp
Fchirp = stft_chirp(m-1, n0-1, M, L, f1, q1);
Fchirp2= Fchirp.^2;   %% square modulus

%% unitary tests:

%% test1: local maximum position
[~,I1] = max(Fchirp2);
[~,I2] = max(spec(:,n0));
if I1 == I2
    fprintf(1,'OK - peak position matches!\n');
else
   fprintf(1,'ERROR - peak position mismatch!\n'); 
end

%% test2: exact amplitude
r = max(Fchirp2) / max(spec(:,n0))
if abs(r-1) < 1e-4
    fprintf(1,'OK - amplitude matches!\n');
else
   fprintf(1,'ERROR - amplitude is wrong: %.2f!\n', r); 
end

plot(m-1, Fchirp.^2, 'r.')   %% display the theoretical value of the window Fh


xlabel('frequency bin')
ylabel('magnitude')
legend('observation', 'model 1 (sinusoid)', 'model 2')  %, 
title(sprintf('M=%d, L=%d, error l2 m1=%.2f, m2=%.2f', M, L, norm(Fh(m-m0,M,L).'-spec(:,n0),2),norm(Fchirp.^2-spec(:,n0),2)))
