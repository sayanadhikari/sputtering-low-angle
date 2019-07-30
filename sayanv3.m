clc; clearvars;  
tar_mat = 'W';
  Observation_mode = input('Enter Observation Mode(e.g. Angle = 0 or Energy = 1) = ');
  switch Observation_mode
      
      case 0
      theta1 =input('Input the first value of incident angle = ');  %angle (in degrees) with normal to target of incident particle.
      theta2 =input('Input the last value of incident angle = '); 
      theta = theta1:1:theta2;
      eo =input(' energy of the incident particle = ');
%       gamma =input('Incident flux of ions (Prt/m^2/s) = ');               
      yldphy = zeros(1,length(theta));
      Ye_E0 = zeros(1,length(theta));
      for i=1:1:length(theta)
%           yldphy(i) = sayandriver(theta(i),eo,z1,z2,am1,am2,es,tgdns);
          yldphy(i)  = sayaninput(theta(i),eo,tar_mat);
      end
  
    %  Converting target density from gms/cm^3 to prt/m^3
    %  1 gram at wt of target has 6.023e23 prt.
    %  => 6.023e26/am2 particles per Kg.
    %  tgdns = tgdns*1.0e3 kg/m^3.
    %  therefore number of prt/m^3 of target = tgdns(kg/m^3)*(6.023e26/am2)(prt/kg)
    
    
    
%       tgdns_ppmc = tgdns*1.0d3*6.023d26/am2;
%       Ytot = yldphy;         % total sputtering yield.
%            Gamma_target = Ytot .* gamma;     %eroded flux from target.
%            Gamma_C = Gamma_target;
%            Erosion_rate = Gamma_C./tgdns_ppmc;   %Erosion rate of target (m/s).

      figure(1);
      plot(theta,real(yldphy));
       grid on;
     xlabel('Theta[degree]');
     ylabel('Physical Sputtering Yield');
     title('Physical Sputtering Yield');
     
     
%      figure(2);
%      loglog(theta,real(Ye_E0));
%      grid on;
%      xlabel('Energy[eV]');
%      ylabel('Physical Sputtering Yield');
%      title('Physical Sputtering Yield Normal incidence');
%      figure(2);
%      figure(2);
%        plot(theta,real(Erosion_rate));
%        grid on;
%      xlabel('Theta[degree]');
%      ylabel('Erosion rate of target');
%      title('Erosion rate of target');
%     %   disp(yldphy);
    
    
    
    case 1
      e1 =input('Input the first value of incident energy = ');  %angle (in degrees) with normal to target of incident particle.
      e2 =input('Input the last value of incident energy = '); 
      ener = linspace(e1,e2,500);
%       ener = [11,12,13,14,15,17,20,25,30,40,50,70,100,140,200,300,500,1000,3000];
%       index1 = find(ener==10);
%       index2 = find(ener==45);
%       index3 = find(ener==100);
%       index4 = find(ener==700);
%       index = [index1 index2 index3 index4];
      theta =input(' incident angle = ');
%       gamma =input('Incident flux of ions (Prt/m^2/s) = ');               
      yldphy = zeros(1,length(ener));
      Ye_E0 = zeros(1,length(ener));
      for i=1:1:length(ener)
          [yldphy(i),Ye_E0(i)] = sayaninput(theta,ener(i),tar_mat);
      end
  
    %  Converting target density from gms/cm^3 to prt/m^3
    %  1 gram at wt of target has 6.023e23 prt.
    %  => 6.023e26/am2 particles per Kg.
    %  tgdns = tgdns*1.0e3 kg/m^3.
    %  therefore number of prt/m^3 of target = tgdns(kg/m^3)*(6.023e26/am2)(prt/kg)
%       tgdns_ppmc = tgdns*1.0d3*6.023d26/am2;
%       Ytot = yldphy;         % total sputtering yield.
%            Gamma_target = Ytot .* gamma;     %eroded flux from target.
%            Gamma_C = Gamma_target;
%            Erosion_rate = Gamma_C./tgdns_ppmc;   %Erosion rate of target (m/s).
     figure(1);
     plot(ener,real(yldphy));
     grid on;
     xlabel('Energy[eV]');
     ylabel('Physical Sputtering Yield');
     title('Physical Sputtering Yield');
     
     figure(2);
     plot(ener,real(Ye_E0));
     grid on;
     xlabel('Energy[eV]');
     ylabel('Physical Sputtering Yield');
     title('Physical Sputtering Yield');
%      figure(2);
%        plot(ener,real(Erosion_rate));
%        grid on;
%      xlabel('Energy[eV]');
%      ylabel('Erosion rate of target');
%      title('Erosion rate of target');
%     %   disp(yldphy);
  end