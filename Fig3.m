clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate Fig.3 of the paper 'Estimation of Instantaneous Frequency and 
% Amplitude of Multi-Component Signals using Sparse Modeling of Signal 
% Innovation', by appling the proposed FRI TLS algorithm onto a speech signal.
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024
%

folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'FRI_lib']));


% Tfr parameter
L = 40;
M = 512;
M2 = floor(M/2);

[s,Fs] = audioread('mData/OBVI.wav'); 

Fs_out = 5000;
s = resample(s,Fs_out,Fs);
N = length(s);

[tfr,stfr]  = tfrsgab2(s, M, L);
Spect = abs(tfr(1:M2,:)).^2;

Spect = Spect(:,300:1000);
Spect(50:end,:)=0;
Spect(1:10,:)=0;

N = size(Spect,2);
t = (0:N-1)/Fs_out;
f = m_axis(M)/M*Fs_out;f = f(1:M/2);



%% Method parameters
Ncomp = 2;                                 % Number of components
M0 = 30;                                   % Frequency truncation - to avoid infinite sum
Method = 2;

Mm = M/2; % frequency support considered for the method (select M/2 if you want to take half of the frequency)
meanF = 125;
m = -meanF:(Mm)-meanF-1;                    % frequency support of the convolution kernel
F = transpose(Fh(m, M, L ));                % Convolution kernel


%% FRI TLS
[tf,~] = estim_FRI(Spect,Ncomp,F,M0,Method,0,0,0);


%% Plot
cols = {'r-.', 'g-.', 'b-.', 'c--', 'm-.', 'g-x', 'w-o'};

Spect = Spect ./max(Spect); % for better spectrogram readability
figure(1)
imagesc(t, f, Spect);
set(gca,'YDir','normal')
colormap gray;
cmap = flipud(colormap);
colormap(cmap);
xlabel('time [s]')
ylabel('frequency [Hz]')
ylim([0 800])

hold on
for c = 1:Ncomp
  IF = tf(:,c)-1;
  h(c) = plot(t,IF/M*Fs_out, cols{c});
  label{c} = sprintf('mode %d', c);
end
legend(h, label);
