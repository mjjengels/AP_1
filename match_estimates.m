function est_sorted = match_estimates(est,trueval)
%   Reorder estimates to match true values.
%
%   Works for small d. Good for this assignment.

est = est(:);
trueval = trueval(:);
d = length(trueval);

P = perms(1:d);

best_err = inf;
best_est = est;

for i = 1:size(P,1)
    candidate = est(P(i,:));
    err = norm(candidate - trueval);
    if err < best_err
        best_err = err;
        best_est = candidate;
    end
end

est_sorted = best_est;
end