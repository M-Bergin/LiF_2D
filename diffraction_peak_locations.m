function [k_out,G,theta_out,phi_out,N_eff]=diffraction_peak_locations(theta_in,phi_in,lambda)

%% Find the fundamental lattice vectors

%k_in=(2*pi)/(5e-11)*[1/sqrt(2),1/sqrt(2)];
%Lattice constant for LiF
a=403.51e-12/sqrt(2);
%a=1;

%Number of peaks to take, better to be odd
N_peaks=15;

n_x=[-floor(N_peaks/2):floor(N_peaks/2)];
G_x=((2*pi)/a)*n_x;

n_y=[-floor(N_peaks/2):floor(N_peaks/2)];
G_y=((2*pi)/a)*n_y;

[G_X,G_Y]=meshgrid(G_x,G_y);
[N_X,N_Y]=meshgrid(n_x,n_y);

G_size=size(G_X);

R=[cosd(phi_in), -sind(phi_in); sind(phi_in), cosd(phi_in)];

Z=R*[G_X(:)';G_Y(:)'];

G_X=reshape(Z(1,:),G_size);
G_Y=reshape(Z(2,:),G_size);

G=sqrt(G_X.^2+G_Y.^2);
N_eff=G*a/(2*pi);%sqrt(N_X.^2+N_Y.^2);

%% Find the outgoing reciprocal vectors

k_mag=(2*pi)/lambda;

k_in=k_mag*[sind(theta_in)*cosd(0),sind(theta_in)*sind(0),cosd(theta_in)];


%Reshape the incoming wavevector
k_in_x=k_in(1)*ones(floor(N_peaks/2)*2+1,1);
k_in_y=k_in(2)*ones(floor(N_peaks/2)*2+1,1);
%k_in_Z=k_in(3)*ones(N_peaks,N_peaks);

[k_in_X,k_in_Y]=meshgrid(k_in_x,k_in_y);

k_out_X=k_in_X+G_X;

k_out_Y=k_in_Y+G_Y;

k_out_Z=sqrt(k_mag^2-(k_out_X.^2+k_out_Y.^2));

k_out=[k_out_X(:),k_out_Y(:),k_out_Z(:)];


%Remove peaks that cannot exist
imag_part=imag(k_out(:,3));
imag_inds=imag_part>0;
k_out(imag_inds,:)=[];
N_eff(imag_inds)=[];

%Find angles
theta_out=acosd(k_out(:,3)/k_mag);
phi_out=atand(k_out(:,2)/k_mag);




