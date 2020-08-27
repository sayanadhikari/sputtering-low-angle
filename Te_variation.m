%% Riemann paper Integration for different alpha
%
%%
clc; clearvars;close all;
%% INITIALIZATION

global alpha c u0 delta

% b = 1;

n = 5000;
multi_alpha = [1 1.25 1.5 1.75 2 2.25 2.5 2.75 3];
me = 0.000548579909; % electron mass in amu
mi = 2.0; % Deuterium mass in amu
% T_ee = [20 40 50 70 100 120 150 200 250];    % Electron temp
% T_ee = linspace(1,200,30);    % Electron temp D-Be
% T_ee = linspace(30,170,30);    % Electron temp D-W  (Temp Ratio 4)
% T_ee = linspace(50,240,5);    % Electron temp D-W  (Temp Ratio 2)
T_ee = linspace(70,300,30);    % Electron temp D-W  (Temp Ratio 1)
temp_ratio = 1;
material = 'W';
% T_ee = 650;
% T_i = T_e/4;    % Ion temp
gamma =1;

% kTebye = T_e;    %kT_e/e
% 
% kT = T_e*1.6E-19;
amu = 2*1.67E-27;
Q = 1.6E-19;
n_CS = 1E19;


%%%%%%%%%%%%%% Progress Bar  %%%%%%%%%%
h = waitbar(0,'1','Name','Phy. Sputtering Calc',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for mm = 1:length(T_ee)
    
    % Check for Cancel button press
    if getappdata(h,'canceling')
        break
    end
    % Report current estimate in the waitbar's message field
    waitbar(mm/length(T_ee),h,sprintf('%2.2f %% done',(100*mm/length(T_ee))))
    
    
    
    T_e = T_ee(mm);
    T_i = T_e*temp_ratio;    % Ion temp
    kTebye = T_e;    %kT_e/e
    kTe = T_e*1.6E-19;
    kTi = T_i*1.6E-19;
    
    M1 = sqrt((2*pi*(me/mi))*(1+(T_i/T_e)));
    for j = 1:length(multi_alpha)
        alpha = multi_alpha(j);
        b(j) = (1/M1)*sind(alpha);
    end
    % a1 = [0.06 0.06 0.06 0.06 0.06];
    % a1 = [0.024 0.03584 0.0475 0.0589 0.0699];
    a1 = [0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06];
    c = 1;

    w = zeros(length(multi_alpha),500);
    II = zeros(length(multi_alpha),500);
    v = zeros(length(multi_alpha),500);
    Chi = zeros(length(multi_alpha),500);
    Vz = zeros(length(multi_alpha),500);
    Vy = zeros(length(multi_alpha),500);
    Vx = zeros(length(multi_alpha),500);
    EE = zeros(length(multi_alpha),500);
    U = zeros(length(multi_alpha),500);
    den = zeros(length(multi_alpha),500);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FOR MARKER %%%%%%%%%%%%%%%
    % markers = {'- *','- o','- s','- v','- +'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k =1:length(multi_alpha)
        
        alpha = multi_alpha(k);
        u0 = c*sind(alpha);
        delta = tand(alpha);
    %     b = (1/M1)*sind(alpha);
        a = linspace(a1(k),b(k),500);
        for j = 1:length(a)
            w(k,j) = ((c+(1/c))/cosd(alpha))-delta*(a(j)+(1/a(j)));
        end
        for i = 1:length(a)
            II(k,i) = -simpsonruleintfun(a(i),b(k),n);
        end
%         fprintf('II= %4.8f \n', II(k,1));
        for l = 1:length(a)
            v(k,l) = sqrt((2*log(a(l)/u0))+(c*c)-(a(l)*a(l))-(w(k,l)*w(k,l)));
        end

        for m = 1:length(a)
            Chi(k,m) = log(a(m));
        end

        Cs = sqrt(gamma*kTi/amu);
%         Cs = sqrt((kTe+(gamma*kTi))/amu);

         for q = 1:length(a)
            Vz(k,q) = Cs*a(q);   
         end

         for r = 1:length(a)
            Vx(k,r) = Cs*w(k,r);
         end

        for s = 1:length(a)
            Vy(k,s) = Cs*v(k,s);
        end

        for n = 1:(length(a)-1)
            EE(k,n) = (Chi(k,(n+1))-Chi(k,(n)))/(II(k,(n+1))-II(k,(n)));
        end

        for o = 1:length(a)
    %         U(o) = -(kTebye)*(log(Vz(o)./Cs));
              U(k,o) = -(kTebye)*(log(a(o)));
        end

        q = 1;
        for p = 1:length(a)
    %         den(p) = (1)*(exp(U(p)/kTebye));
            den(k,p) = n_CS.*(1)/(a(p));
        end
        den_norm = den(k,:)./max(den(k,:));
    % den = den./0.4348E21;
    %     den = den_norm.*n_CS;

        Flux(k) = den(k,length(a))*Vz(k,length(a));


        theta(k) = atand(Vz(k,length(a))/sqrt((Vx(k,length(a))*Vx(k,length(a)))+(Vy(k,length(a))*Vy(k,length(a)))));


%         KE(k) = (1/(1.6E-19))* (0.5*amu*( (Vx(k,length(a)).^2)+ (Vy(k,length(a)).^2) + (Vz(k,length(a)).^2)));
% %         KE(k) = KE(k)./3;
%         KE3d(mm,k) = KE(k);

             E0(k) = (0.5*(Cs^2)*amu/Q)-(log(Cs*sind(alpha))*(amu/Q));
             KEx(k) = Vx(k,length(a))^2 ;
             KEy(k) = Vy(k,length(a))^2 ;
             KEz(k) = Vz(k,length(a))^2 ;
             CONST = ((0.5*amu)/Q); 
             KE(k) =  CONST*(KEx(k) + KEy(k) + KEz(k));
             fprintf('KEx =  %e ,KEy =  %e, KEz =  %e, KE =  %e, CONST = %e, E0 = %e\n',KEx(k),KEy(k),KEz(k),KE(k),CONST,E0);
             KE3d(mm,k) = KE(k);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%    SPUTTERING PART 1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Phys(k) = sayaninput(theta(k),KE(k),material);

        Phys3d(mm,k) = Phys(k);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    end
    % h(1).Marker = '*';
    % h(2).Marker = 'o';
    % h(3).Marker = 's';
    % h(4).Marker = '+';
    % h(5).Marker = 'v';

%         hold all;
%         figure(1);
%     %     subplot(2,2,1);
%     %     plot(II,a,markers{mod(k,numel(markers))+1});
% %         plot(II(1,:),a,'-',II(2,:),a,'-',II(3,:),a,'-',II(4,:),a,'-',II(5,:),a,'-');
%         plot(II(5,:),a,'-');
%         xlabel(' -z/\rho_i');
%         ylabel('u');
%         legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');


    %     hold all;
    %     figure(2);
    % %     subplot(2,2,2);
    %     plot(II,a,'-');
    %     xlabel(' -z/\rho_i');
    %     ylabel('u');
    %     legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');
    %     
    %     hold all;
    %     figure(3);
    % %     subplot(2,2,3);
    %     plot(II,v,'-');
    %     xlabel(' -z/\rho_i');
    %     ylabel('v');
    %     legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');
    %     
    % %     hold all;
    % %     
    % %     subplot(2,2,4);
    % %     plot(II,Chi,'-');
    % %     xlabel(' -z/\rho_i');
    % %     ylabel('\chi');
    %     %legend('\alpha = 1','\alpha = 1.5','\alpha = 2','\alpha = 2.5','\alpha = 3');
    %     hold all;
    %     figure(4);
    % %     subplot(2,2,1);
    %     plot(II(1:(length(II)-1)),abs(EE),'-');
    %     xlabel(' -z/\rho_i');
    %     ylabel('E');
    %     legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');
    %     
    %     hold all;
    %     figure(5);
    % %     subplot(2,2,2);
    %     plot(II,U,'-');
    %     xlabel(' -z/\rho_i');
    %     ylabel('U');
    %     legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');
    %     
    %     hold all;
    %     figure(6);
    % %     subplot(2,2,3);
    %     plot(II,den_norm,'-');
    %     xlabel(' -z/\rho_i');
    %     ylabel('n');
    %     legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');

    % legend('\alpha = 1^{\circ}','\alpha = 1.5^{\circ}','\alpha = 2^{\circ}','\alpha = 2.5^{\circ}','\alpha = 3^{\circ}');
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%    SPUTTERING PART 2     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
%         tgdns = 1.85; % For Be
%         am2 = 9.012;    % For Be

%      
      am2 = 183.8; % For W
      tgdns = 19.3; %For W
%         
        
          tgdns_ppmc = tgdns*1.0d3*6.023d26/am2;
          Ytot = Phys;         % total sputtering yield.
          for nn=1:length(multi_alpha)
              Gamma_target(nn) = Ytot(nn) .* Flux(nn);     %eroded flux from target.
              Gamma_C = Gamma_target;
              Erosion_rate(nn) = Gamma_C(nn)./tgdns_ppmc;   %Erosion rate of target (m/s).
          end
    % 
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % %         figure(1);
%         hold all;
%         figure(11);
%     %     subplot(2,2,1);
%         plot(multi_alpha,Flux,'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
%         xlabel('Grazing Angle(\alpha^{\circ})');
%         ylabel('Flux(m^{-1}.sec^{-1})','Interpreter', 'tex');
%         grid on;
%         legend('T_i = T_e','T_i = T_e/2','T_i = T_e/4');
%         figure(7);
%     %     subplot(2,2,2);
%         plot(multi_alpha,theta,'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
%         xlabel('Grazing Angle(\alpha^{\circ})');
%         ylabel('Incident Angle(\theta^{\circ})');
%         grid on;
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %     subplot(2,2,3);
%         figure(5)
%         hold all;
%         plot(multi_alpha,KE,'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
%         xlabel('Grazing Angle(\alpha^{\circ})');
%         ylabel('Incident Ion Energy(eV)');
%         grid on;
%         legend('T_i = T_e','T_i = T_e/2','T_i = T_e/4');
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         hold all;
%         figure(2);
%     %     subplot(2,2,4);
%         plot(multi_alpha,Phys,'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
%         xlabel('Angle(\alpha^{\circ})');
%         ylabel('Physical Sputtering Yield');
%         grid on;
%         legend('T_e = 10','T_e = 30','T_e = 45','T_e = 50');
    %     
    %     
%         hold all;
%         figure(3);
%         plot(multi_alpha,Erosion_rate,'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
%         xlabel('Angle(\alpha^{\circ})');
%         ylabel('Erosion Rate(m.sec^{-1})');
%     %     title('Erosion Rate');
%         grid on;
%         legend('T_e = 10','T_e = 30','T_e = 45','T_e = 50');
        
        

end
delete(h)       % DELETE the waitbar; don't try to CLOSE it.


figure(1)
% surf(multi_alpha,T_ee,Phys3d,'EdgeColor','none');
surf(multi_alpha,T_ee,Phys3d,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud');
% ylim([20 250]);
% ylim([20 700]);
view(70,26);
% xlabel('Angle(\alpha^{\circ})','rot',-60,'pos',[-57,2416,-0.032]);
% ylabel('Electron Temperature(eV)','rot',15,'pos',[3.545,77.531,-0]);
xlabel('Angle(\alpha^{\circ})','rot',-45);
ylabel('Electron Temperature(eV)','rot',7);
title('Physical Sputtering Yield');
% camlight; lighting phong;


% figure(6)
% plot(T_ee,KE3d(:,1),'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
% %         ylim([0 80]);
% %         view(84,26);
% xlabel('Electron Temperature(eV)');
% ylabel('Ion Impact Energy(eV)');
% title('Variation of Ion Impact Energy');
% grid on;

% 
% figure(7)
% plot(T_ee,KE3d(:,1),'-o','LineWidth',1.5,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','r',...
%                     'MarkerSize',5);
% %         ylim([0 80]);
% %         view(84,26);
% xlabel('Electron Temperature(eV)');
% ylabel('Ion Impact Energy(eV)');
% title('Variation of Ion Impact Energy');
% grid on;


% figure(6)
% surf(multi_alpha,T_ee,KE3d,'EdgeColor','none');
% ylim([0 80]);
% view(84,26);
% xlabel('Angle(\alpha^{\circ})');
% ylabel('Electron Temperature(eV)');
% title('Ion Impact Energy(eV)');
