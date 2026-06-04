function [w,t,theta_grid,response] = minimax_beamformer(M,Delta,theta1,theta2,theta_grid)

if nargin < 5
    theta_grid = -90:0.1:90;
end

Agrid = steervec(M, Delta, theta_grid);
a1 = steervec(M, Delta, theta1);
a2 = steervec(M, Delta, theta2);

cvx_begin quiet
    variable w(M) complex
    variable t nonnegative

    minimize(t)

    subject to
        w' * a1 == 1
        w' * a2 == 0
        abs(w' * Agrid) <= t
cvx_end

response = abs(w' * Agrid);

end