%% Calculation of fitting values f,b,c for D-Be (ion-target)
function [fy,by,cy,theta0star]= fitting_parameter_D_Be(eo)

E0 = [11,12,13,14,15,17,20,25,30,40,50,70,100,140,200,300,500,1000,3000];
f = [13.5302,13.1593,13.3940,12.5003,11.7442,11.2812,10.6565,10.1966,10.0263,9.3179,9.1102,8.4246,7.6933,7.3337,6.6082,6.3383,6.1678,5.8465,5.8554];
b = [10.2782,9.6932,9.5880,8.7410,7.9688,7.2662,6.4178,5.7049,5.3137,4.5386,4.2137,3.5596,2.9730,2.6468,2.1695,1.9598,1.7951,1.5024,1.2982];
c = [0.7530,0.4965,0.5074,0.5374,0.5699,0.6167,0.6964,0.7740,0.8192,0.8830,0.9064,0.9434,0.9653,0.9660,0.9733,0.9468,0.8947,0.8107,0.7851];
Esp = 1.0;
theta0star = 180-acosd(sqrt(1/(1+(eo/Esp))));
% theta0star = theta0star*1.7453293d-2;

fit_data_f = fit(E0',f','smoothingspline');
fit_data_b = fit(E0',b','smoothingspline');
fit_data_c = fit(E0',c','smoothingspline');


fy = feval(fit_data_f,eo);
by = feval(fit_data_b,eo);
cy = feval(fit_data_c,eo);
end