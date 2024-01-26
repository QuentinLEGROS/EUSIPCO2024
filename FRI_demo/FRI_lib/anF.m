function h = anF(Method,y,NbDiracs,M,boolnoise)


if Method == 1 % Yulle-Walker system
    h = YW(y);
elseif Method == 2 % TLSM
    h = AF_TLSM(y,NbDiracs,M,boolnoise)';
else
    error('Select Method in [1,2]');
end
