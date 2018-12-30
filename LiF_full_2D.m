% Function to approximate the diffraction pattern from a LiF crystal and
% what signal it would produce in a SHeM

%Parameters
lambda=6.63e-34/(sqrt(5*4*1.67e-27*1.38e-23*(273+25)));
theta_in=45;
phi_in=45-59;
alph=0.5; %Flux in diffuse component

%Detector aperture dimensions 
x_ap_mid=sqrt(2)*3e-3;
y_ap_mid=0;
r_ap=0.5e-3;%0.5e-3;

%Calculate positions of the diffraction peaks
[k_out,G,theta_out,phi_out,N_eff]=diffraction_peak_locations(theta_in,phi_in,lambda);

%Generate k space grid to plot diffraction pattern on
k_x_max=1.1e11;
k_x_N=1000;

k_y_max=1.1e11;
k_y_N=1000;

k_x_vec=linspace(-k_x_max,k_x_max,k_x_N);
k_y_vec=linspace(-k_y_max,k_y_max,k_y_N);

[k_X,k_Y]=meshgrid(k_x_vec,k_y_vec);

%Set k_Z by energy conservation
k_mag=(2*pi)/lambda;
k_Z=sqrt(k_mag^2-(k_X.^2+k_Y.^2));
imag_inds=find(imag(k_Z)>0);

I_k_unrot=zeros(k_x_N,k_y_N);

N_points=size(k_out,1);

%Set the width of the peaks and how they decay
peak_width=3.5e9;
pattern_width=2;

%Add in the diffraction pattern
for n_point=1:N_points
    I_x_temp=normpdf(k_x_vec,k_out(n_point,1),peak_width);
    I_y_temp=normpdf(k_y_vec,k_out(n_point,2),peak_width);
    if N_eff(n_point)==0
        I_k_unrot=I_k_unrot+(normpdf(N_eff(n_point),0,pattern_width)*(I_y_temp'*I_x_temp));
    elseif abs(N_eff(n_point)-1)<0.5
        I_k_unrot=I_k_unrot+(normpdf(N_eff(n_point),0,pattern_width)*(I_y_temp'*I_x_temp));    
    else
        I_k_unrot=I_k_unrot+(normpdf(N_eff(n_point),0,pattern_width)*(I_y_temp'*I_x_temp));
    end
end

%Normalise the diffraction pattern contribution
I_k_unrot=(I_k_unrot/(sum(sum(I_k_unrot))))*(1-alph);


%Add in diffuse component
theta_mat=real(atand(sqrt(k_X.^2+k_Y.^2)./k_Z));
I_diff=real(cosd(theta_mat));

%Normalise
I_diff=(I_diff./(sum(sum(I_diff))))*alph;

%Calculate total diffraction pattern with diffuse component.
I_tot=I_k_unrot+I_diff;
%I_tot=I_diff;

%Remove the imaginary parts
I_tot(imag_inds)=0;


I_k=I_tot;
I_diff_rot=I_diff;


%figure;imagesc(k_y_vec,k_x_vec,I_k)



%Find the k values that are entering the aperture


z_vec=linspace(1e-3,8e-3,100); %3e-3/sqrt(2);
N_z=length(z_vec);

I_z=NaN*zeros(N_z,1);
I_z_diff=NaN*zeros(N_z,1);

%Set up plotting
figure;
subplot(1,2,1)
imagesc((I_k))
set(gca,'Yticklabel',[])
set(gca,'Xticklabel',[])
set(gca,'YTick',[]);
set(gca,'XTick',[]);
axis equal tight
xlabel('k_x')
ylabel('k_y')
hold on
h_p1=plot(0,0,'r.');
subplot(1,2,2)
h_p2=plot(z_vec,I_z);
hold on
h_p3=plot(z_vec,I_z);
xlabel('Sample distance/m')
ylabel('Relative intensity')
xlim([z_vec(1) z_vec(end)])
%ylim([0 0.2])

for n_z=1:N_z
    
    z=z_vec(n_z);
    
    %Find the position on the plate if the crystal wasn't rotated
    X_pos_unrot=real(k_X./k_Z)*z;
    Y_pos_unrot=real(k_Y./k_Z)*z;
    
    %Shift if using unrotated coordinates
    X_pos_unrot=X_pos_unrot+z;
    
    %Clear the outsides that shouldn't be included
    imag_inds=find(imag(k_Z)>0);
    X_pos_unrot(imag_inds)=NaN;
    Y_pos_unrot(imag_inds)=NaN;
    
    %Rotate the co-ordinates
    X_pos_rot=X_pos_unrot*cosd(-phi_in)-Y_pos_unrot*sind(-phi_in)+z;
    Y_pos_rot=X_pos_unrot*sind(-phi_in)+Y_pos_unrot*cosd(-phi_in);
    
    
    
    %Find the distance from the aperture for the k-space matrix
    ap_dist_unrot=sqrt(((X_pos_unrot-x_ap_mid)/sqrt(2)).^2+((Y_pos_unrot-y_ap_mid)).^2);
    
    %Rotate this matrix to match the actual co-ordinates
    %ap_dist=imrotate(ap_dist_unrot,phi_in);
    ap_dist=ap_dist_unrot;
    
    nan_ind=find(ap_dist==0);
    ap_dist(nan_ind)=NaN;
    
    %Clip the output pattern for only the transmitted flux
    ind_in=find(ap_dist<r_ap);
    
    I_clip=NaN*I_k;
    I_clip_diff=NaN*I_k;
    I_bound=NaN*I_k;
    I_clip(ind_in)=I_k(ind_in);
    I_clip_diff(ind_in)=I_diff_rot(ind_in);
    I_bound(ind_in)=1;
    
    I_b=imbinarize(I_bound);
    I_edge=edge(I_b);
    [edge_ind1,edge_ind2]=find(I_edge==1);
    
    
    %bw=bwboundaries(I_b);
    %bw_plot=bw{1};
    
    %hold on
    %plot(bw_plot(:,2),bw_plot(:,1),'r','LineWidth',1)
    %plot(edge_ind2,edge_ind1,'rx')
    
    h_p1.XData=edge_ind2;
    h_p1.YData=edge_ind1;
    
    I_z(n_z)=nansum(nansum(I_clip));
    I_z_diff(n_z)=nansum(nansum(I_clip_diff));
    
    h_p2.YData=I_z;
    h_p3.YData=I_z_diff;
    drawnow
end


%theta=atand((3e-3*sqrt(2)-z_vec)./z_vec);

% % Create thesis figures
% 
% figure;imagesc(k_y_vec,k_x_vec,I_k)
% xlabel('k_x /m^{-1}')
% ylabel('k_y /m^{-1}')
% axis equal tight
% set(gca,'FontSize',12,'LineWidth',1)
% 
% 
% 
% z=3e-3/sqrt(2);
% 
% %Find the position on the plate if the crystal wasn't rotated
% X_pos_unrot=real(k_X./k_Z)*z;
% Y_pos_unrot=real(k_Y./k_Z)*z;
% 
% %Shift if using unrotated coordinates
% X_pos_unrot=X_pos_unrot+z;
% 
% %Clear the outsides that shouldn't be included
% imag_inds=find(imag(k_Z)>0);
% X_pos_unrot(imag_inds)=NaN;
% Y_pos_unrot(imag_inds)=NaN;
% 
% %Rotate the co-ordinates
% X_pos_rot=X_pos_unrot*cosd(-phi_in)-Y_pos_unrot*sind(-phi_in)+z;
% Y_pos_rot=X_pos_unrot*sind(-phi_in)+Y_pos_unrot*cosd(-phi_in);
% 
% 
% 
% %Find the distance from the aperture for the k-space matrix
% ap_dist_unrot=sqrt(((X_pos_unrot-x_ap_mid)/sqrt(2)).^2+((Y_pos_unrot-y_ap_mid)).^2);
% 
% %Rotate this matrix to match the actual co-ordinates
% %ap_dist=imrotate(ap_dist_unrot,phi_in);
% ap_dist=ap_dist_unrot;
% 
% nan_ind=find(ap_dist==0);
% ap_dist(nan_ind)=NaN;
% 
% %Clip the output pattern for only the transmitted flux
% ind_in=find(ap_dist<r_ap);
% 
% I_clip=NaN*I_k;
% I_clip_diff=NaN*I_k;
% I_bound=NaN*I_k;
% I_clip(ind_in)=I_k(ind_in);
% I_clip_diff(ind_in)=I_diff_rot(ind_in);
% I_bound(ind_in)=1;
% 
% I_b=imbinarize(I_bound);
% I_edge=edge(I_b);
% [edge_ind1,edge_ind2]=find(I_edge==1);
% 
% 
% hold on
% 
% N_I=size(I_k);
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
% % print(['C:\Users\mberg\Dropbox\Thesis\Contrast\Figures\LiF_diffract_pattern.pdf'],'-dpdf','-r0')
% % savefig('C:\Users\mberg\Dropbox\Thesis\Contrast\Figures\LiF_diffract_pattern.fig')
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
% xlim([1 8])
% 
% 
% legend('Full signal','Diffuse signal')
% 
% % print -depsc2 C:\Users\mberg\Dropbox\Thesis\Contrast\Figures\LiF_simul_z.eps
% % savefig('C:\Users\mberg\Dropbox\Thesis\Contrast\Figures\LiF_simul_z.fig')
% 




