function A = steervec(M, Delta, theta_deg)
%   A = steervec(M, Delta, theta_deg)
%
%   M         : number of antennas
%   Delta     : antenna spacing in wavelengths
%   theta_deg : d x 1 or 1 x d vector of DOAs in degrees
%
%   A is M x d, with columns a(theta_i).

theta_deg = theta_deg(:).';        % row vector
m = (0:M-1).';                     % antenna indices

A = exp(1j * 2*pi*Delta * m * sind(theta_deg));
end