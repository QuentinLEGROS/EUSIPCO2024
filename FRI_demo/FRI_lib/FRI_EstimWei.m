function [FoundDiracsWeights] = FRI_EstimWei(mv,FoundDiracsLocations,NbDiracs,LenF,y,modF,m)
% [FoundDiracsWeights] = FRI_EstimWei(mv,FoundDiracsLocations,NbDiracs,LenF,y)
%
% Weight retrivial by solving a linear system
%
% INPUT:
% mv                        : Indices of Fourier coefficient kept
% FoundDiracsLocations      : Dirac positions
% NbDiracs                  : Number of Dirac to retrieve
% LenF                      : Signal length
% y
%
% OUTPUT:
% FoundDiracsWeights        : Dirac Weights
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 12-may-2021

pas = 1;
p=0;
for hh = -pas:pas
    p=p+1;
    FoundDiracsLocations = sort(mod(FoundDiracsLocations - (LenF-modF) -1, LenF))-1+pas;

    C = kron(mv',round(FoundDiracsLocations(1:NbDiracs)));
    C = exp((-1i * 2 * pi * C) ./ LenF);
    weights(:,p) = sqrt(2*abs(lsqr(C,y)));
end
FoundDiracsWeights = mean(weights,2)

% FoundDiracsWeights = sqrt(2*abs(lsqr(C,y)));

