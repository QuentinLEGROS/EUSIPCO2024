function [ z ] = Gk(omega, omega0,T, k)

z = 1 ./ (1 + 1j*(omega-omega0)*T).^k;

end

