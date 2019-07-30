%% Calculation of fitting values f,b,c for D-W (ion-target)
function [fy,by,cy,theta0star]= fitting_parameter_D_W(eo)

E0 = [250,270,300,350,400,500,600,700,1000];
f = [3.9269,3.3257,2.4382,2.2277,1.9710,1.5306,1.0222,1.2293,1.2531];
b = [2.5577,1.9540,1.4079,1.2410,1.0176,0.6843,0.3585,0.3613,0.2141];
c = [0.8242,0.9277,0.9691,0.9674,0.9898,1.0177,1.0601,1.0416,1.0543];
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