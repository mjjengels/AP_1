function [w,t,theta_grid,response] = minimax_beamformer(M,Delta,theta1,theta2,theta_grid)
%   Computes a minimax beamformer for one desired direction.
%   The beamformer is constrained to have unit response at theta1:
%
%       w^H a(theta1) = 1
%
%   It is also constrained to place a null at theta2:
%
%       w^H a(theta2) = 0
%
%   The maximum response over theta_grid is minimized using CVX.

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
