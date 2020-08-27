 
% c  Inputs to the subroutine
% c  theta -> angle (in degrees) with normal to target of incident particle.
% c  eo	-> energy of the incident particle.
% c  z1	-> nuclear charge of incident atom.
% c  z2	-> nuclear charge of target atom.
% c  am1	-> atomic mass of incident atom (amu).
% c  am2	-> atomic mass of target atom (amu).
% c  es	-> surface binding energy (heat of sublimation) of target (eV).
% c  tgdns -> target density (gms/cc)
% 
% c  Output from the subroutine
% c  yldphy -> The physical sputtering yield.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yldphy,Ye_E0] = sputtInput(theta,eo,tar_mat)
        e0=1.6E-19; %electron charge
%   theta =input('Input the incident angle = ');  %angle (in degrees) with normal to target of incident particle.
%   eo =input(' energy of the incident particle = ');
%   z1 =input('nuclear charge of incident atom = ');
%   z2=input('nuclear charge of target atom = ');
%   am1 =input('atomic mass of incident atom (amu) = ');
%   am2 =input('atomic mass of target atom (amu) = ');
%   es =input(' surface binding energy (heat of sublimation) of target (eV) = ');
%   tgdns =input('target density (gms/cc) = ');
% c  Initialising
%%%%%% Extra for Material Choosing %%%%%
    % c  Inputs to the subroutine
%   Incident_atom = input('Enter Incident Atom = ');
%   switch Incident_atom
%       case 'H'
%           z1 = 1;
%           am1 = 1.008;
%       case 'D'
          z1 = 1;
          am1 = 2.016;
%       otherwise
%         %warning('Unexpected Incident atom not in current database');
%         exit;
%   end
   Target_atom = tar_mat;
  switch Target_atom
      case 'W'
          z2 = 74;
          am2 = 183.8;
          es = 8.9;
          tgdns = 19.3;     %19.3;
          % Fitting parameters for D-W (ion-target)
          lambda = 3.665e+0;
          qtotal = 1.0e-4;
          mu = 2.279;
          Eth =  229.743;   %es * ( (7.0 * (am2/am1)^(-0.54)) + (0.15 * (am2/am1)^(1.12)));            %9.090;
          nu = 1.210;
          [fy,by,cy,theta0star]= fitting_parameter_D_W(eo);
      case 'Be'
          z2 = 4;
          am2 = 9.012;
          es = 3.32;
          tgdns = 1.85;  %   1850; %kg/m^3     %1.85 g/cc
          % Fitting parameters for D-Be (ion-target)
          lambda = 7.991e+0;
          qtotal = 2.204e-3;
          mu = 2.946;
          Eth =  9.090;   %es * ( (7.0 * (am2/am1)^(-0.54)) + (0.15 * (am2/am1)^(1.12)));            %9.090;
          nu = 1.072;
          [fy,by,cy,theta0star]= fitting_parameter_D_Be(eo);
      otherwise
        %warning('Unexpected Target atom not in current database');
         exit;
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     eps_L = 2.82110e+2;
%     lambda = 1.7575;
%     qtotal = 0.1044;
%     mu = 1.9906;
%     Eth = 9.5059;
%     nu = 1.0;
%     eps_L  = 2.82110e+2;
 %%%%%%%%%%% calculation of epsilon_L %%%%%%%%%%%%%%%%%%%%
    pwr1by3 = 1/3;
    pwr2by3 = 2/3;
    pi_term = ((9*pi^2)/128)^pwr1by3;
    a_B = 0.0529177;
    z123z223sum = (z1^pwr2by3)+(z2^pwr2by3);
    a_L = pi_term*a_B*(1/sqrt(z123z223sum));
    eps_L = eo*(am2/(am1+am2)) * a_L/(z1*z2*e0^2);    %2.82110e+2;
    CONST_term = pi_term*a_B/e0^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       etf = 30.74d0 .* (am1+am2)./am2  .* z1.*z2.*z123z223sum.^0.5d0; %%%30.74d0
       

      eobyetf = eo ./ etf;
      stoppwr1 = (0.5d0.*log(1.0d0+1.2288d0.*eobyetf));
      stoppwr2 = (eobyetf + 0.1728d0.*eobyetf.^0.5d0 + 0.008d0.*eobyetf.^0.1504d0);
      stoppwr = stoppwr1 ./ stoppwr2;
    
    
    n = 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yldphy,Ye_E0] = sputtYeilds(theta,eo,lambda,qtotal,mu,Eth,stoppwr,eobyetf,nu,eps_L,n,fy,by,cy,theta0star,z1,z2,am1,am2,es);
% disp(yldphy);
end
