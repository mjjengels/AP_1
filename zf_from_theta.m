function [W,Ahat] = zf_from_theta(M,Delta,theta_hat)
%   W satisfies approximately W^H Ahat = I.

Ahat = steervec(M, Delta, theta_hat);
W = Ahat / (Ahat' * Ahat);
end