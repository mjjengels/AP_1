function [W,Ahat,Shat,theta_hat,t,theta_grid,Y] = minimax_from_freq(X,f_hat,Delta,theta_grid)

if nargin < 3 || isempty(Delta)
    Delta = 1/2;
end

if nargin < 4 || isempty(theta_grid)
    theta_grid = -90:0.1:90;
end

[M,N] = size(X);
f_hat = f_hat(:);
d = length(f_hat);

if d > M
    error('Need d <= M for zero-forcing.');
end

k = 0:N-1;
Shat = exp(1j * 2*pi * f_hat * k);

A_ls = X * pinv(Shat);

m = (0:M-1).';
theta_hat = zeros(d,1);

for i = 1:d
    ai = A_ls(:,i);

    if abs(ai(1)) > 0
        ai = ai / ai(1);
    end

    ph = unwrap(angle(ai));
    p = polyfit(m, ph, 1);

    sin_theta = p(1) / (2*pi*Delta);
    sin_theta = max(min(real(sin_theta),1),-1);

    theta_hat(i) = asind(sin_theta);
end

Ahat = steervec(M, Delta, theta_hat);

theta_grid = unique([theta_grid(:); theta_hat(:)]).';
Agrid = steervec(M, Delta, theta_grid);

W = zeros(M,d);
t = zeros(d,1);
Y = zeros(d,length(theta_grid));

for i = 1:d
    c = zeros(d,1);
    c(i) = 1;

    cvx_begin quiet
        variable w(M) complex
        variable ti nonnegative

        minimize(ti)

        subject to
            Ahat' * w == c
            abs(Agrid' * w) <= ti
    cvx_end

    W(:,i) = w;
    t(i) = ti;
    Y(i,:) = abs(w' * Agrid);
end

end