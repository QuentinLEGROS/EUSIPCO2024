function [FoundDiracsLocations] = FRI_EstimLoc(h,NbDiracs,LenF,modF)
% [FoundDiracsLocations] = FRI_EstimLoc(h,NbDiracs,LenF,modF)
%
% Position retrivial through the computation of the roots of the Z
% transform of h
%
% INPUT:
% h             : annihilating filter
% NbDiracs      : Number of Dirac to retrieve
% LenF          : Signal length
% modF          : Mode of the postulated data distribution
%
% OUTPUT:
% FoundDiracsLocations    : Position(s) estimée(s)
%
% Author: Q.Legros (quentin.legros@telecom-paris.fr)
% Date: 12-may-2021

% Method: annihilating filter search root of the z-transform (closest to unit circle)
TZroot = conj(roots(h));
[~,index] = sort(abs(abs(TZroot)-1));
TZroot = TZroot(index(1:NbDiracs));
% Location

% FoundDiracsLocations = sort(mod( LenF *  imag(log(TZroot.'))  ./ (2*pi), LenF));
% FoundDiracsLocations = FoundDiracsLocations - (LenF-modF) -1;

FoundDiracsLocations = LenF * imag(log(TZroot.'))  ./ (2*pi);



