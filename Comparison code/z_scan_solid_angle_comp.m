% Function to compare the signal obtained from a diffuse z scan to the
% solid angle of the detector aperture

%Load HOPG data for comparison. 
load('Z000013.mat')

%Zero position in mm
z_zero=(5195e3-0.9e6)/1e6; %HOPG
%z_zero=3.575; % \pm 0.15-0.2mm Flat
%z_zero=(5195e3-0.9e6)/1e6;

%Experimental z positions in mm and m
z2=z_zero-z_position/1e6;
z3=z2*1e-3;

%% Plot z scan
fig_h=figure;
yyaxis left
plot(z2,current_z,'x','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',15,'LineWidth',1)
xlabel('z/mm')
ylabel('Measured current/A')

yyaxis right



%Parameter setup for solid angle calculation
h_max=8e-3;
h_N=1000;
h_delta=h_max/h_N;

h_vec=h_delta:h_delta:h_max;

Omega=NaN*zeros(h_N,1);


for n=1:h_N
    
    %Set the parameters for the ellipse at this z position
    h=h_vec(n);
    a=sqrt(2)*0.5e-3;
    b=0.5e-3;
    p=3e-3*sqrt(2)-h;
    q=0;
    
   %Perform the analytic integral
    integrand= @(phi) (1- h./sqrt(p.^2+q.^2+h.^2+2*a*p*cos(phi)+2*b*q*sin(phi)+...
        a^2*cos(phi).^2+b^2*sin(phi).^2)).*((a*b+p*b*cos(phi)+q*a*sin(phi))./...
        (p^2+q^2+2*a*p*cos(phi)+2*b*q*sin(phi)+a^2*cos(phi).^2+b^2*sin(phi).^2));
    
    Omega(n)=integral(integrand,0,2*pi);
end

%Plot the result
plot(h_vec/1e-3,Omega,'LineWidth',1.5)
hold on
plot([3/sqrt(2) 3/sqrt(2)],[0 0.1-1e-7],'k--','LineWidth',1.5)
xlabel('z/mm')
ylabel('\Omega')

set(gca,'FontSize',15,'LineWidth',1)

xlim([1 8])
%  print('C:\Users\mberg\Dropbox\LiF\z_scan_solid_angle_comp.eps','-depsc2')
%  savefig('C:\Users\mberg\Dropbox\LiF\z_scan_solid_angle_comp.fig')

