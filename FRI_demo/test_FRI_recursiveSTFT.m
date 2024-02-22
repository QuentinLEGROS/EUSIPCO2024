clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Unitary test FRI on recusive STFT
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer (dominique@fourer.fr)
%  Date    : 24-feb-2024
%


folder = './';
%% required paths 
addpath(folder);
addpath(strcat([folder 'tools']));
addpath(strcat([folder 'synchrosqueezedSTFT']));
addpath(strcat([folder 'FRI_lib']));
addpath(strcat([folder 'RecursiveSTFT']));
addpath(strcat([folder 'modulation']));
addpath(strcat([folder 'Compute_Amplitude_DF']));
%% Time-frequency representation parameters
N     = 500;       %% signal length
M     = 500;       %% Number of frequential bin
L     = 20;        %% analysis window size (in bin)
npad  = 50;



%% Define signal x0
Amp(:,1) = 1*ones(N,1);


X(:,1) = Amp(:,1).*(fmconst(N, 0.1));
tf0(:,1) = 0.1*ones(N,1);

Ncomp = size(X,2);                          %% number of components
x0 = sum(X,2);
X = transpose(X);



%% Method parameters
M0 = 10;                                    % Frequency truncation - to avoid infinite sum
Method = 2;
% tf = zeros(M/2,N); ia = tf;                 % Init
tf = zeros(N,Ncomp); ia = tf;                 % Init
k=3;                                        % recursive filter order


%% Compute reference
[tfr]  = tfrgab2(x0, M, L);
Spect0=(abs(tfr(1:round(M/2),:))).^2;


%% Noise
SNR = 20;
x = sigmerge(x0, randn(size(x0)), SNR);



%% Define analysis window and parameters for the recursif STFT
[F,a,b] = init_recursif_data(M,L,k);


xp = [zeros(k-1,1);x];
tfrp    = zeros(M/2, N+k);   % Initialization

figure(1)
f_range = (0:(floor(M/2)-1))/M;
n_range = 1:N;
%imagesc(n_range, f_range, Spect);
set(gca,'YDir','normal')
colormap('gray');cmap = colormap;cmap = flipud(cmap);
colormap(cmap);
for n = k:N+k-1
    tfrp(:,n+1) = transpose(sum(b.*xp(n-k+1:n),1)) - sum(a.*tfrp(:,n-k+1:n),2);

    Spect(:,n-k+1) = abs(tfrp(:,n+1)).^2;
    % Spect(:,n-k+1) = abs(tfrp(:,n+1));
    
if mean(Spect(:,n-k+1))>=1e-3
    [tf(n-k+1,:)] = estim_FRI_recursif(Spect(:,n-k+1),Ncomp,F,M0,Method);
else
    tf(n-k+1,:) = 1;
    ia(n-k+1,:) = 1;
end
imagesc(n_range(1:(n-k+1)), f_range, Spect(:,1:(n-k+1)));
hold on
% imagesc(abs(tfrp).^2)
plot(tf(1:(n-k+1))/M,'k')
hold off
xlabel('Time index');
ylabel('Normalized Frequency');
pause(0.01)
end



%% Estimation at nt = 250
nt = 250;
figure(2);
plot(f_range, Spect(:,nt));
hold on;
plot(tf(nt)/M, Spect(round(tf(nt))-1,nt), 'rx');
legend(sprintf('Spectrogram at n=%d',nt), sprintf('estimated IF (m/M = %.3f)',tf(nt)/M));

%% TODO amplitude estimation

%% Le reste ne sert
% 
% [tfr] = tfrgab2(x, M, L);
% Amp = Oracle_Amp_DF(tfr,Ncomp, L,round(tf));
% %% plot Window - spectrogram
% Spect = abs(tfrp(:,k+1:end)).^2;
% 
% figure(1)
% hold on
% % plot(Spect(:,100)./sum(Spect(:,100)),'r')
% % plot(circshift(F,-77)./sum(circshift(F,-77)),'k')
% plot(Spect(:,100),'r')
% plot(circshift(F,-77),'k')
% hold off
% 
% % 
% %% Main
% % [tfr0] = tfrgab2(x, M, L); %% compute SST
% % Spect0 = abs(tfr0(1:M/2,:)).^2;
% % [tf,ia] = estim_FRI(Spect,Ncomp,F,M0,Method,L,tf0(:,1).*M);
% 
% 
% 
% % Plots
% figure(1)
% subplot(2,1,1)
% imagesc((Spect0))
% subplot(2,1,2)
% plot(tf)
% ylim([0 250])
% 
% figure(2)
% title('amplitudes')
% % subplot(2,1,1)
% hold on
% plot(ia(npad:end-npad,1))
% plot(Amp(npad:end-npad,1))
% hold off
% legend('estimation','Ground truth')
% % subplot(2,1,2)
% % hold on
% % plot(ia(npad:end-npad,2))
% % plot(Amp(npad:end-npad,2))
% % hold off
% % legend('estimation','Ground truth')

