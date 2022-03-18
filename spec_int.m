function [output] = spec_int(mu, v)
%SPEC_INT specitla integral
%   mu is scalar, v is vector
output = 1/(2*mu) - v/mu*sqrt(pi/mu).*exp(v.^2/mu).*qfunc(v*sqrt(2/mu));
end

