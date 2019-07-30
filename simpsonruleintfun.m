%% SIMPSON 1/3 RULE FOR NUMERICAL INTEGRATION
%
function I = simpsonruleintfun(a,b,n)

%% INITIALIZATION
h = (b-a)/n;

sum1 = 0.0;
for i = 1:2:(n-1)
    xi = a+(i*h);
    sum1 = sum1 + intfunc(xi);
end

sum2 = 0.0;
for i = 2:2:(n-1)
    xi = a+(i*h);
    sum2 = sum2 + intfunc(xi);
end

I = (h/3)*(intfunc(a)+intfunc(b)+(4*sum1)+(2*sum2));
% fprintf('I= %4.8f \n', I)
end
    