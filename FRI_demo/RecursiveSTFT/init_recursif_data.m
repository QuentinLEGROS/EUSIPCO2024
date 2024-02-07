function [F,a,b] = init_recursif_data(M,L,k)


%% module
omega0 = 2*pi * (M/4)/M;
omega = -(M/4)+1:(M/4);
T = L;%L/(M/2);
% F = (1+((omega-omega0)*T).^2).^(-k/2);
F = (1+((omega-omega0)*T).^2).^(-k);
F = F./max(F);


%% Recursive_STFT
mm = 1:M/2;
lambda = (1+mm-1)/M;
pTs    = -1.0/L + 1i*2*pi.*lambda;
alpha  = exp(pTs);
[a,b]  = Gk22(k, L, alpha);
b = flip(transpose(b),1); a = flip(a(:,2:end),2);



