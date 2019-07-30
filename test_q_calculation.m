
clc; clearvars;
z1 = 1;
am1 = 2.016;
z2 = 4;
am2 = 9.012;
es = 3.38;
tgdns = 1.85;
e0=1.6E-19; %electron charge
eo = 500;

pwr1by6 = 1.0/6.0;
pwr5by6 = 5.0/6.0;
pwr1by3 = 1.0/3.0;
pwr2by3 = 2.0/3.0;

z123xz223 = (z1^pwr2by3)*(z2^pwr2by3);
z123z223sum = (z1^pwr2by3)+(z2^pwr2by3);
am2byam1 = (am2/am1);
% Expression for q       
qpar1 = (1.633 * z123xz223 * z123z223sum^pwr1by3 )/(es^pwr2by3);
qpar2 = (am1^pwr5by6 * am2^pwr1by6) / (am1 + am2);
qpar3 = (0.15 + (0.05*am2byam1)) /(1.0 + (0.05*am2byam1^1.6));
qtotal = (qpar1 * qpar2 * qpar3);
 pi_term = ((9*pi^2)/128)^pwr1by3;
a_B = 5.29e-11;
a_L = pi_term*a_B*(1/sqrt(z123z223sum));
eps_L = eo*(am2/(am1+am2)) * a_L/(z1*z2*e0^2); 
CONST_term = e0*pi_term*a_B/e0^2;

