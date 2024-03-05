function [FoundDiracsLocations] = FRI_EstimLoc(h,NbDiracs,LenF)
%
% Position retrivial through the computation of the roots of the Z
% transform of h
%
% INPUT:
% h             : annihilating filter
% NbDiracs      : Number of Dirac to retrieve
% LenF          : Signal length
%
% OUTPUT:
% FoundDiracsLocations    : estimated IF
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024


TZroot = conj(roots(h));
[~,index] = sort(abs(abs(TZroot)-1));
TZroot = TZroot(index(1:NbDiracs));

% Location
FoundDiracsLocations = LenF * imag(log(TZroot.'))  ./ (2*pi);
