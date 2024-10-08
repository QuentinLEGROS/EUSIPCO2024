function [FoundDiracsWeights] = FRI_EstimWei(mv,FoundDiracsLocations,NbDiracs,LenF,y)
% [FoundDiracsWeights] = FRI_EstimWei(mv,FoundDiracsLocations,NbDiracs,LenF,y)
%
% Weight retrivial by solving a linear system
%
% INPUT:
% mv                        : Indices of Fourier coefficient kept
% FoundDiracsLocations      : Dirac positions
% NbDiracs                  : Number of Dirac to retrieve
% LenF                      : Signal length
%
% OUTPUT:
% FoundDiracsWeights        : estimated IA
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024



% FoundDiracsLocations = round(FoundDiracsLocations);
C = kron(mv',FoundDiracsLocations(1:NbDiracs)); 
C = exp((-1i * 2 * pi * C) ./ LenF);

FoundDiracsWeights = sqrt(abs(lsqr(C,y)))*norm(y);
end





