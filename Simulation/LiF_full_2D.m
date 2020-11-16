% Function to approximate the diffraction pattern from a LiF crystal and
% what signal it would produce in a SHeM

%%%%%%%%%%%% Setup of parameters %%%%%%%%%%%

%Parameters
lambda=6.63e-34/(sqrt(5*4*1.67e-27*1.38e-23*(273+25)));
theta_in=45;
phi_in=0;%45-59;
alph=0.5; %Flux in diffuse component

%Detector aperture dimensions
x_ap_mid=3.5e-3*2;%sqrt(2)*3e-3;
y_ap_mid=0;
r_ap=0.25e-3;

%Offset of aperture from pinhole plate
z_offset=1.5e-3;
%z positions to calculate
z_min=1.5e-3; %z in m
z_max=12e-3;
N_z=200; %Number of points, increase to 1000 for better quality


%%%%%%%%%% Creation of scattering distribution in k space %%%%%%%%%%

%Calculate positions of the diffraction peaks
[k_out,G,theta_out,phi_out,N_eff]=diffraction_peak_locations(theta_in,phi_in,lambda);

%Generate k space grid to plot diffraction pattern on
k_x_max=1.2e11;
k_x_N=1000;

k_y_max=1.2e11;
k_y_N=1000;

k_x_vec=linspace(-k_x_max,k_x_max,k_x_N);
k_y_vec=linspace(-k_y_max,k_y_max,k_y_N);

[k_X,k_Y]=meshgrid(k_x_vec,k_y_vec);

%Set k_Z by energy conservation
k_mag=(2*pi)/lambda;
k_Z=sqrt(k_mag^2-(k_X.^2+k_Y.^2));
imag_inds=imag(k_Z)>0;



%Create the diffraction pattern
I_k=zeros(k_x_N,k_y_N);

N_points=size(k_out,1);

%Set the width of the peaks and how they decay
peak_width=3.5e9;
pattern_width=2;

%Main loop to add in the diffraction pattern
for n_point=1:N_points
    I_x_temp=normpdf(k_x_vec,k_out(n_point,1),peak_width);
    I_y_temp=normpdf(k_y_vec,k_out(n_point,2),peak_width);
    if N_eff(n_point)==0
        I_k=I_k+(normpdf(N_eff(n_point),0,pattern_width)*(I_y_temp'*I_x_temp));
    elseif abs(N_eff(n_point)-1)<0.5
        I_k=I_k+(normpdf(N_eff(n_point),0,pattern_width)*(I_y_temp'*I_x_temp));
    else
        I_k=I_k+(normpdf(N_eff(n_point),0,pattern_width)*(I_y_temp'*I_x_temp));
    end
end

%Normalise the diffraction pattern contribution
I_k=(I_k/(sum(sum(I_k))))*(1-alph);


%Add in diffuse component
theta_mat=(atand(sqrt(k_X.^2+k_Y.^2)./k_Z));
theta_mat(theta_mat~=real(theta_mat))=NaN;

%Distribution for diffuse scattering
I_diff=~isnan(theta_mat);
%Distribution for the solid angle
%I_diff=1./real(cosd(theta_mat));

%Normalise the diffuse contribution
I_diff=(I_diff./(nansum(nansum(I_diff))))*alph;


%Calculate total diffraction pattern with diffuse component.
I_tot=I_k+I_diff;


%Remove the imaginary parts
I_tot(imag_inds)=0;



%%%%%%%% Integration of distribution to find signal %%%%%%%%%%


%Setup
z_vec=linspace(z_min,z_max,N_z);
I_z=NaN*zeros(N_z,1);
I_z_diff=NaN*zeros(N_z,1);
k_corr_fact=NaN*zeros(N_z,1);

%Set up plotting
f_h=figure;
subplot(1,2,1)
imagesc((I_tot))
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[])
set(gca,'YTick',[]);
set(gca,'XTick',[]);
axis equal tight
xlabel('k_x')
ylabel('k_y')
hold on
h_p1=plot(0,0,'r.','MarkerSize',12);
%set(gca,'FontSize',24,'LineWidth',3)
subplot(1,2,2)
h_p2=plot(z_vec-z_offset,I_z,'LineWidth',3);
hold on
h_p3=plot(z_vec-z_offset,I_z,'LineWidth',3);
xlabel('z/m')
ylabel('Relative intensity')
xlim([z_vec(1) z_vec(end)]-z_offset)
ylim([0 0.08])
%set(gca,'FontSize',24,'LineWidth',3)
%f_h.WindowState='fullscreen';

%Main loop
for n_z=1:N_z
    
    %Set the z position
    z=z_vec(n_z);
    
    %Find the position on the plate for each k value
    X_pos=real(k_X./k_Z+1)*z;
    Y_pos=real(k_Y./k_Z)*z;
    
    
    %Clear the outsides that shouldn't be included
    imag_inds=find(imag(k_Z)>0);
    X_pos(imag_inds)=NaN;
    Y_pos(imag_inds)=NaN;
    
    
    %Find the distance from the aperture for the k-space matrix
    ap_dist=sqrt(((X_pos-x_ap_mid)/sqrt(2)).^2+((Y_pos-y_ap_mid)).^2);
    %figure;imagesc(abs(ap_dist));caxis([0 0.05])
    
    %Remove
    nan_ind=(ap_dist==0);
    ap_dist(nan_ind)=NaN;
    
    %Clip the output pattern for only the transmitted flux
    ind_in=find(ap_dist<r_ap);
    
    I_clip=NaN*I_tot;
    I_clip_diff=NaN*I_tot;
    I_bound=NaN*I_tot;
    I_clip(ind_in)=I_tot(ind_in);
    I_clip_diff(ind_in)=I_diff(ind_in);
    I_bound(ind_in)=1;
    
    
    %Find the edge of the detection aperture
    I_b=imbinarize(I_bound);
    I_edge=edge(I_b);
    [edge_ind1,edge_ind2]=find(I_edge==1);
    
    %Plot the new aperture position
    h_p1.XData=edge_ind2;
    h_p1.YData=edge_ind1;
    
    %Integrate the collected flux in the aperture
    I_z(n_z)=nansum(nansum(I_clip));
    I_z_diff(n_z)=nansum(nansum(I_clip_diff));
    
    %Debug plot
    %figure;imagesc(temp2)
    
    %Update plot
    h_p2.YData=I_z;
    h_p3.YData=I_z_diff;
    drawnow
    
    %Create a video
    %F(n_z) = getframe(gcf);
    
    
    
    %bw=bwboundaries(I_b);
    %bw_plot=bw{1};
    %hold on
    %plot(bw_plot(:,2),bw_plot(:,1),'r','LineWidth',1)
    %plot(edge_ind2,edge_ind1,'rx')
end

%theta=atand((3e-3*sqrt(2)-z_vec)./z_vec);


% write the frames to file
% for i=1:2:length(F)
%     imwrite(F(i).cdata,['Temp/',num2str(i),'.png'])
% end




%%%%%%%%%%%% Create thesis figures %%%%%%%%%%%%%%%%%%%%
%Ignore if you aren't re-making my thesis.

% figure;imagesc(k_y_vec,k_x_vec,I_k)
% xlabel('k_x /m^{-1}')
% ylabel('k_y /m^{-1}')
% axis equal tight
% set(gca,'FontSize',12,'LineWidth',1)


%Repeat above process at specular only

z=x_ap_mid/2;%3e-3/sqrt(2);

%Find the position on the plate for each k value
X_pos=real(k_X./k_Z+1)*z;
Y_pos=real(k_Y./k_Z)*z;


%Clear the outsides that shouldn't be included
imag_inds=find(imag(k_Z)>0);
X_pos(imag_inds)=NaN;
Y_pos(imag_inds)=NaN;


%Find the distance from the aperture for the k-space matrix
ap_dist=sqrt(((X_pos-x_ap_mid)/sqrt(2)).^2+((Y_pos-y_ap_mid)).^2);
%figure;imagesc(abs(ap_dist));caxis([0 0.05])


nan_ind=find(ap_dist==0);
ap_dist(nan_ind)=NaN;

%Clip the output pattern for only the transmitted flux
ind_in=find(ap_dist<r_ap);

I_clip=NaN*I_tot;
I_clip_diff=NaN*I_tot;
I_bound=NaN*I_tot;
I_clip(ind_in)=I_tot(ind_in);
I_clip_diff(ind_in)=I_diff(ind_in);
I_bound(ind_in)=1;


%Plot the edge of the detection aperture
I_b=imbinarize(I_bound);
I_edge=edge(I_b);
[edge_ind1,edge_ind2]=find(I_edge==1);


N_I=size(I_k);

% hold on
% plot(edge_ind2*(max(k_x_vec)-min(k_x_vec))/N_I(1)+min(k_x_vec),edge_ind1*(max(k_y_vec)-min(k_y_vec))/N_I(2)+min(k_y_vec),'r.')
%
% fig_h=gcf;
%
%
% set(fig_h,'Units','Inches');
% pos = get(fig_h,'Position');
% set(fig_h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% set(gca,'FontSize',15,'LineWidth',1)
%
% % print(['C:\Users\mberg\Dropbox\LiF\LiF_diffract_pattern.pdf'],'-dpdf','-r0')
% % savefig('C:\Users\mberg\Dropbox\LiF\LiF_diffract_pattern.fig')
%
%
% % z scan plot
%
% figure;plot(z_vec/1e-3,I_z,'LineWidth',1)
% hold on
% plot(z_vec/1e-3,I_z_diff,'LineWidth',1)
%
% xlabel('Sample distance/mm')
% ylabel('Relative intensity')
%
% set(gca,'FontSize',15,'LineWidth',1)
% xlim([z_min z_max]/1e-3)
%
%
% legend('Full signal','Diffuse signal')
%
% % print -depsc2 C:\Users\mberg\Dropbox\LiF\LiF_simul_z.eps
% % savefig('C:\Users\mberg\Dropbox\LiF\LiF_simul_z.fig')


%%%%%%%%%%%%%%%% Adjust for new pinhole plate %%%%%%%


z_vec=z_vec-z_offset;

z_min=z_min-z_offset;
z_max=z_max-z_offset;


%%%%%%%%%%%%%%%% Create paper figure %%%%%%%%%%%%%%%%%

% plot on large axes
figure;plot(z_vec/1e-3,I_z,'LineWidth',1.5)
hold on
plot(z_vec/1e-3,I_z_diff,'LineWidth',1.5)
%plot([z z]/1e-3,[0 max(I_z)],'k--','LineWidth',1.5)

xlabel('z/mm')
ylabel('Relative intensity')

set(gca,'FontSize',15,'LineWidth',1)
xlim([z_min z_max]/1e-3)

% create smaller axes in top right, and plot on it
axes('Position',[.53 .5 .4 .4])
box on

imagesc(k_y_vec/1e11,k_x_vec/1e11,I_tot)
xlh =xlabel('k_x /10^{11}m^{-1}');
xlh.Position=[0.1 1.6 1];
ylh=ylabel('k_y /10^{11}m^{-1}');
ylh.Position=[-1.75 0 1];
axis equal tight
set(gca,'FontSize',15,'LineWidth',1)
set(gca,'XTick',[-1,-0.5,0,0.5,1]);
set(gca,'YTick',[-1,-0.5,0,0.5,1]);

hold on

plot((edge_ind2*(max(k_x_vec)-min(k_x_vec))/N_I(1)+min(k_x_vec))/1e11,(edge_ind1*(max(k_y_vec)-min(k_y_vec))/N_I(2)+min(k_y_vec))/1e11,'r.','MarkerSize',4)

p1 = [0.4658,0.534];%[0.3579 0.1046];                         % First Point
%p2 = [0.3028 -0.1046];                         % Second Point
p2 = [0.2037,-0.5494];%[0.2522 -0.3248];                         % Second Point
dp = p2-p1;                         % Difference

%quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',2,'MaxHeadSize',1,'Color','r')
%text(0.12,-0.43,'[011]','FontSize',14,'Color','r');
%text(-0.07,-0.75,'[011]','FontSize',15,'Color','r');

fig_h=gcf;


set(fig_h,'Units','Inches');
pos = get(fig_h,'Position');
set(fig_h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gca,'FontSize',14,'LineWidth',1)


% print(['C:\Users\mberg\Dropbox\LiF\LiF_diffract_pattern_new.pdf'],'-dpdf','-r1000')
% savefig('C:\Users\mberg\Dropbox\LiF\LiF_diffract_pattern_new.fig')

 %print('C:\Users\mab679\Dropbox\Newcastle\LiF_mid_3_5.png','-dpng','-r1000')


