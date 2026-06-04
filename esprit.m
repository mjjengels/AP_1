function theta = esprit(X,d)
%   Estimate DOAs from X using spatial ESPRIT.
%   theta = esprit(X,d)
%
%   X : M x N data matrix
%   d : number of sources
%
%   theta is returned in degrees.

[M,~] = size(X);

if d >= M
    error('For spatial ESPRIT, need d < M.');
end

% Signal subspace from SVD
[U,~,~] = svd(X, 'econ');
Us = U(:,1:d);

% Selection matrices for two shifted subarrays
J1 = [eye(M-1), zeros(M-1,1)];
J2 = [zeros(M-1,1), eye(M-1)];

% Shift-invariance equation:
% J2*Us = J1*Us*Psi
Psi = pinv(J1*Us) * (J2*Us);

lambda = eig(Psi);

% For ULA:
% lambda_i = exp(j*2*pi*Delta*sin(theta_i))
% Assignment uses Delta = 1/2.
Delta = 1/2;
sin_theta = angle(lambda) / (2*pi*Delta);

% Numerical clipping
sin_theta = max(min(real(sin_theta),1),-1);

theta = asind(sin_theta);
theta = sort(theta(:));
end