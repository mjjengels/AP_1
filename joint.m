function [theta,f] = joint(X,d,m)
%   X : M x N received data
%   d : number of sources
%   m : time-smoothing factor
%
%   Constructs the extended matrix
%
%       Xe = [X(:,1:L);
%             X(:,2:L+1);
%             ...
%             X(:,m:N)]
%
%   where L = N-m+1.
%
%   This gives a Khatri-Rao structured model:
%
%       Xe = (B_f o A_theta) S0
%
%   with spatial and temporal shift invariances.

[M,N] = size(X);

if m < 2
    error('m must be at least 2.');
end

L = N - m + 1;

if L < d
    error('Need N-m+1 >= d.');
end

if d >= M
    error('Need d < M for spatial shift invariance.');
end

% Build time-smoothed extended data matrix
Xe = zeros(M*m, L);

for ell = 1:m
    rows = (ell-1)*M + (1:M);
    Xe(rows,:) = X(:, ell:ell+L-1);
end

% Signal subspace of extended matrix
[Ue,~,~] = svd(Xe, 'econ');
E = Ue(:,1:d);

% Spatial shift selection inside each time block
Jm_theta_1 = kron(eye(m), [eye(M-1), zeros(M-1,1)]);
Jm_theta_2 = kron(eye(m), [zeros(M-1,1), eye(M-1)]);

Psi_theta = pinv(Jm_theta_1*E) * (Jm_theta_2*E);

% Temporal shift selection between consecutive smoothed blocks
Jm_f_1 = kron([eye(m-1), zeros(m-1,1)], eye(M));
Jm_f_2 = kron([zeros(m-1,1), eye(m-1)], eye(M));

Psi_f = pinv(Jm_f_1*E) * (Jm_f_2*E);

Psi(:,:,1) = Psi_theta;
Psi(:,:,2) = Psi_f;

[D,~] = joint_diag(Psi);

lambda_theta = diag(D(:,:,1));
lambda_f = diag(D(:,:,2));

Delta = 1/2;

sin_theta = angle(lambda_theta) / (2*pi*Delta);
sin_theta = max(min(real(sin_theta),1),-1);

theta = asind(sin_theta);
f = mod(angle(lambda_f)/(2*pi), 1);

% Sort by theta while preserving theta-frequency pairing
[theta,idx] = sort(real(theta(:)));
f = real(f(idx));
end


function [D,T] = joint_diag(Psi)

[d,~,K] = size(Psi);

C = zeros(d,d);

for k = 1:K
    nk = norm(Psi(:,:,k), 'fro');
    if nk > 0
        C = C + (1j)^(k-1) * Psi(:,:,k) / nk;
    end
end

[T,~] = eig(C);

D = zeros(d,d,K);

for k = 1:K
    D(:,:,k) = T \ Psi(:,:,k) * T;
end

end