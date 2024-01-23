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
% y
%
% OUTPUT:
% FoundDiracsWeights    : Dirac Weights
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 12-may-2021

C = kron(mv',FoundDiracsLocations(1:NbDiracs));
C = exp((-1i * 2 * pi * C) ./ LenF);

FoundDiracsWeights = abs(transpose(pinv(C)*y));