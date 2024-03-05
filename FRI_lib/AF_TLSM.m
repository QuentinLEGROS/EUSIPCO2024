function H = AF_TLSM(y,NbDiracs)
% H = AF_TLSM(y,NbDiracs,M)
%
% Algorithme Annihilating filter: total least-squares method présenté dans
% l'article 'Sparse Sampling of Signal Innovations: Theory, Algorithms and Performance Bounds'
%
% INPUT:
% y : data
% NbDiracs : nombre de Diracs
% M : troncature utilisé pour cadzow denoising
%
% OUTPUT:
% H : filtre annihilateur
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 1-mar-2024

L = NbDiracs;
%% Toeplitz matrix
for j=0:L
    T(:,j+1)=y(L-j+1:end-j);
end

% SVD
[~,~,V] = svd(T);

% Estimation of H : vecteur propre correspondant a la plus petite valeur
% propre de S
H = conj(V(:,end));

    




