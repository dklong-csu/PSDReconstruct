function y = Ispline(x,k,i,t)
%Mspline TODO
%   TODO
arguments
    x {mustBeNumeric,mustBeFinite};
    k {mustBeInteger,mustBePositive};
    i {mustBeInteger,mustBePositive};
    t {mustBeFloat};
end

y=0*x;
for iii=1:numel(x)
    y(iii) = integral(@(y) Mspline(y,k,i,t),t(1),x(iii));
end



end