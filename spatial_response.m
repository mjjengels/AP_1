function y = spatial_response(w,M,Delta,theta_grid)
%SPATIAL_RESPONSE Compute y(theta) = |w^H a(theta)|.

Agrid = steervec(M, Delta, theta_grid);
y = abs(w' * Agrid);
end