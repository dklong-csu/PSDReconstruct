function y = Mspline(x,k,i,t)
%Mspline TODO
%   TODO
arguments
    x {mustBeNumeric,mustBeFinite};
    k {mustBeInteger,mustBePositive};
    i {mustBeInteger,mustBePositive};
    t {mustBeFloat};
end
if (t(i+k) - t(i)) < 1e-8
    y = 0*x;
elseif k == 1
    y = (x >= t(i)) .* (x < t(i+1)) * 1/( t(i+1) - t(i) );
else
    coeff = k/ ( (k-1) * ( t(i+k)-t(i) ) );
    recursiveM1 = (x - t(i)) .* Mspline(x,k-1,i,t);
    recursiveM2 = (t(i+k) - x) .* Mspline(x,k-1,i+1,t);
    y = coeff * (recursiveM1 + recursiveM2);
end



end