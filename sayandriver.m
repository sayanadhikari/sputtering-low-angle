% c  Subroutine to evaluate sputtering yeilds of
% c  any target material due to Physical Sputtering.
% 
% c  NOTE : These formula are valid only for normal incidence
% c         and have to be modified for angular incidence.
% 
% c  REF-1: "Erosion Processes in plasma-wall Interactions"
% c         C. Garci'a Rosales
% c	  Jnl. Nucl. Mat. 211 (1994) 202-214

% subroutine physput(theta,eo,z1,z2,am1,am2,es,tgdns,yldphy)
% implicit real*8(a-h,o-z)

% real*8 eo, z1, z2, am1, am2, es, eth, etf, E_L, stoppwr, yeildpar, ethbyeo, yldphy
 
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
function [yldphy,Ye_E0] = sayandriver(theta,eo,lambda,qtotal,mu,Eth,stoppwr,eobyetf,nu,eps_L,n,fy,by,cy,theta0star,z1,z2,am1,am2,es)
% c  Inputs to the subroutine

%   z1 =input('nuclear charge of incident atom = ');
%   z2=input('nuclear charge of target atom = ');
%   am1 =input('atomic mass of incident atom (amu) = ');
%   am2 =input('atomic mass of target atom (amu) = ');
%   es =input(' surface binding energy (heat of sublimation) of target (eV) = ');
%   tgdns =input('target density (gms/cc) = ');
% c  Initialising

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
      qtotal1 = (qpar1 * qpar2 * qpar3);

      
% Expression for E0/Eth

      eobyeth = eo/Eth;%eo/Eth;

% Nuclear stopping cross-section calculation based on Kr-C potential
    
    w_eps_L = eps_L + (0.1728*sqrt(eps_L)) + (0.008 * eps_L^0.1504);
    
    Sn_Kr_C = stoppwr; %( (0.5 * log(1+(1.2288*eps_L)) )/( (w_eps_L)^nu) )^n;

%Sputtering Tield for normal incident

    eobyeth_minus_one_pwr_mu =((eobyeth-1)^mu);
    
    Ye_E0 = qtotal*Sn_Kr_C*(eobyeth_minus_one_pwr_mu/(lambda+eobyeth_minus_one_pwr_mu));
    

% Angular Contrubution term
%     theta = theta*1.7453293d-2;
    thetabythetaostar = theta/theta0star;
    
    cosine_term = (cos((thetabythetaostar*90*1.7453293d-2)^cy));
    
    theta_term = (((cosine_term)^-fy)*exp(by*(1-(1/cosine_term))));
    
    Ye_E0_theta = Ye_E0*theta_term;


% c  Calculation of the physical sputtering yeild
      yldphy = Ye_E0_theta;
      if (yldphy<=0.0d0) 
          yldphy = 0.0d0;
      end
end