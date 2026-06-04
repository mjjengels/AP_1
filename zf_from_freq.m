function [W,Ahat,Shat] = zf_from_freq(X,f_hat)
%   First reconstruct Shat from estimated frequencies.
%   Then estimate Ahat by least squares:
%
%       Ahat = X Shat^H (Shat Shat^H)^(-1)
%
%   Finally compute W = Ahat (Ahat^H Ahat)^(-1).

[M,N] = size(X);

k = 0:N-1;
Shat = exp(1j * 2*pi * f_hat * k);

Ahat = X * Shat' / (Shat * Shat');

W = Ahat / (Ahat' * Ahat);
end