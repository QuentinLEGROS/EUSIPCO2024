function h = anF(Method,y,NbDiracs)

% Apply either the Prony method of alternative to compute / estimate the
% annihiliating filter
%
% INPUT:
% Method     : For replacing the Prony method in the presence of noise (1:without | 2:TLSD)
% y          : filter to annihilate
% NbDiracs   : Number of components
% M          : frequency support
%
%
%
% OUTPUT:
% h           : annihilating filter
%
%  Authors : Q.Legros (quentin.legros@univ-orleans.fr) and D.Fourer
%  Date    : 1-mar-2024



if Method == 1 % Yulle-Walker system
    h = YW(y);
elseif Method == 2 % TLSM
    h = AF_TLSM(y,NbDiracs)';
else
    error('Select Method in [1,2]');
end
