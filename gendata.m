function [X,A,S] = gendata(M,N,Delta,theta,f,SNR)
%   X     : M x N received data
%   A     : M x d array response matrix
%   S     : d x N source matrix
%   M     : number of antennas
%   N     : number of samples
%   Delta : antenna spacing in wavelengths
%   theta : d x 1 source directions in degrees
%   f     : d x 1 normalized source frequencies
%   SNR   : SNR per source in dB

theta = theta(:);
f = f(:);
d = length(theta);

if length(f) ~= d
    error('theta and f must have the same length.');
end

% Array response
A = steervec(M, Delta, theta);

% Source matrix
k = 0:N-1;
S = exp(1j * 2*pi * f * k);    % d x N

% Noiseless data
Xclean = A*S;

% Each source has unit power. SNR per source = 1 / noise_variance.
SNRlin = 10^(SNR/10);
sigma2 = 1/SNRlin;
sigma = sqrt(sigma2);

Noise = sigma/sqrt(2) * (randn(M,N) + 1j*randn(M,N));

X = Xclean + Noise;
end