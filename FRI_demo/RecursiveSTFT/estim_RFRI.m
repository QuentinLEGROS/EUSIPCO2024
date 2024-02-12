function [tf] = estim_RFRI(x,Fr,M,N,k,a,b,Ncomp,Method,M0)




xp = [zeros(k-1,1);x];
tfrp    = zeros(M/2, N+k);   % Initialization
% Spect = zeros(N,M);

for n = k:N+k-1
    tfrp(:,n+1) = transpose(sum(b.*xp(n-k+1:n),1)) - sum(a.*tfrp(:,n-k+1:n),2);
    Spect(:,n-k+1) = abs(tfrp(:,n+1)).^2;
    if sum(Spect(:,n-k+1))>=1e-6
        [tf(n-k+1,:)] = estim_FRI_recursif(Spect(:,n-k+1),Ncomp,Fr,M0,Method);
    else
        tf(n-k+1,:) = NaN;
        % ia(n-k+1,:) = NaN;
    end
end











