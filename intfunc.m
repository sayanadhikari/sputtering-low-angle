function y = intfunc(u)
global alpha c u0 delta
f = (c*c)+(2*log(u/u0))-(u*u)-(((c+(1/c))/cosd(alpha))-(delta*(u+(1/u))))^2;
y = ((1-(u*u))/u)*(1/sqrt(f));
end