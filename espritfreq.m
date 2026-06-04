function f = espritfreq(X,d)
%   Estimate normalized frequencies from X.
%   f = espritfreq(X,d)
%
%   Uses the time Vandermonde structure in S.
%   Since SVD gives the right singular vectors V, and V spans
%   the conjugated temporal subspace, the eigenvalues are approximately
%   exp(-j*2*pi*f_i).

[~,N] = size(X);

if d >= N
    error('For frequency ESPRIT, need d < N.');
end

% Right signal subspace
[~,~,V] = svd(X, 'econ');
Vs = V(:,1:d);

% Time-shift selection
J1 = [eye(N-1), zeros(N-1,1)];
J2 = [zeros(N-1,1), eye(N-1)];

Psi = pinv(J1*Vs) * (J2*Vs);

lambda = eig(Psi);

% Because Vs corresponds to S^H:
% lambda_i ≈ exp(-j*2*pi*f_i)
f = mod(-angle(lambda)/(2*pi), 1);

f = sort(real(f(:)));
end