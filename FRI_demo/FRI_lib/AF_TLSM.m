function H = AF_TLSM(y,NbDiracs,M,boolnoise)
% H = AF_TLSM(y,NbDiracs,M,boolnoise)
%
% Algorithme Annihilating filter: total least-squares method présenté dans
% l'article 'Sparse Sampling of Signal Innovations: Theory, Algorithms and Performance Bounds'
%
% INPUT:
% y : data
% NbDiracs : nombre de Diracs
% M : troncature utilisé pour cadzow denoising
% boolnoise : booleen sur la quantité de bruit : application ou non de
% Cadzow denoising
%
% OUTPUT:
% H : filtre annihilateur
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 12-may-2021

% booleen pour application de Cadzow denoising
if boolnoise==1
%     y = Cadzow_denoising(y,NbDiracs,M);
    %% Rectangular toeplitz matrix
%     for j=0:M
%         T(:,j+1)=y(M-j+1:end-j);
%     end
%     
%     T = Cadzow_denoising_upgraded(T,NbDiracs);
%     
    V = Cadzow_denoising(y,NbDiracs,M);
    H = conj(V(:,end));
else
    L = NbDiracs;
    %% Toeplitz matrix
    for j=0:L
        T(:,j+1)=y(L-j+1:end-j);
    end
    % SVD
    [U,S,V] = svd(T);
    % Estimation of H : vecteur propre correspondant a la plus petite valeur
    % propre de S
    H = conj(V(:,end));
end
    




