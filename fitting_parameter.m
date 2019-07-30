%% Calculation of fitting values f,b,c for D-Be (ion-target)
clc; clearvars;
eo = 33; %incident energy

tar_mat = 'W';

Target_atom = tar_mat;
switch Target_atom
      case 'Be'
      E0 = [11,12,13,14,15,17,20,25,30,40,50,70,100,140,200,300,500,1000,3000];
      f = [13.5302,13.1593,13.3940,12.5003,11.7442,11.2812,10.6565,10.1966,10.0263,9.3179,9.1102,8.4246,7.6933,7.3337,6.6082,6.3383,6.1678,5.8465,5.8554];
      b = [10.2782,9.6932,9.5880,8.7410,7.9688,7.2662,6.4178,5.7049,5.3137,4.5386,4.2137,3.5596,2.9730,2.6468,2.1695,1.9598,1.7951,1.5024,1.2982];
      c = [0.7530,0.4965,0.5074,0.5374,0.5699,0.6167,0.6964,0.7740,0.8192,0.8830,0.9064,0.9434,0.9653,0.9660,0.9733,0.9468,0.8947,0.8107,0.7851];
      
      case 'W'
      E0 = [250,270,300,350,400,500,600,700,1000];
      f = [3.9269,3.3257,2.4382,2.2277,1.9710,1.5306,1.0222,1.2293,1.2531];
      b = [2.5577,1.9540,1.4079,1.2410,1.0176,0.6843,0.3585,0.3613,0.2141];
      c = [0.8242,0.9277,0.9691,0.9674,0.9898,1.0177,1.0601,1.0416,1.0543];
      otherwise
        %warning('Unexpected Target atom not in current database');
         exit;
 end
Esp = 1.0;
theta0 = 180-acosd(sqrt(1/(1+(eo/Esp))));
% figure(1);
% % hold on;
% plot(E0,f,'o -- b');
fit_data_f = fit(E0',f','smoothingspline');
fit_data_b = fit(E0',b','smoothingspline');
fit_data_c = fit(E0',c','smoothingspline');


figure(1);
plot(E0,f);
xlim( [0, 1000] )
hold on;
plot(fit_data_f);
fy = feval(fit_data_f,eo);
plot(eo,fy,'*','LineWidth',2)
hold off;

figure(2);
plot(E0,b);
xlim( [0, 1000] )
hold on;
plot(fit_data_b);
by = feval(fit_data_b,eo);
plot(eo,by,'*','LineWidth',2)
hold off;

figure(3);
plot(E0,c);
xlim( [0, 1000] )
hold on;
plot(fit_data_c);
cy = feval(fit_data_c,eo);
plot(eo,cy,'*','LineWidth',2)
hold off;

% figure(1);
% plot(fit_data_f,E0,f);
% hold on;
% fy = feval(fit_data_f,eo);
% plot(eo,fy,'*','LineWidth',2)
% 
% figure(2);
% plot(fit_data_b,E0,b);
% hold on;
% by = feval(fit_data_b,eo);
% plot(eo,by,'*','LineWidth',2)
% 
% figure(3);
% plot(fit_data_c,E0,c);
% hold on;
% cy = feval(fit_data_c,eo);
% plot(eo,cy,'*','LineWidth',2)