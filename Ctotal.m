function Ctotal()

%load excel crystal plasticity finite element data
A = xlsread('G120x120x10_Ngrains594_random_tensionX_inc300.xlsx','B:D');
B = xlsread('G120x120x10_Ngrains594_random_tensionX_inc300','F:H');
C = [A, B];

%create axis and hydrostatic stress arrays[Pa]
x=C(:,1); 
y=C(:,2);
z=C(:,3); 
delta0=((C(:,4)+C(:,5)+C(:,6))/3.0)*10^6; 
length(x)

%number of 3D meshes 90*100*10
Nx=max(x)+1; 
Ny=max(y)+1;
Nz=max(z)+1;
hs=zeros(Nx,Ny,Nz);

%create a 3D matrix of hydrostatic stress
for i=1:size(C,1) 
    hs(x(i)+1,y(i)+1,z(i)+1)=delta0(i);
end
eval([ 'vtkwrite(''hydrostastic_stress.vtk'',' '''structured_points'','  '''hydrostastic_stress'', hs) ' ]) 

format long;


% calculate steps [m]
dx=1.0e-6;
dy=1.0e-6;
dz=1.0e-6;

%material property parameters
nprint=200; %total output
D= 2.88e-10; %hydrogen atom diffusion parameter [m2 s-1]
v_H = 1.67e-6; %patitial molar mass of hydrogen atom [m3 mol-1]
R = 8.314; %gas constant [N m kmol-1K-1] [m3 ¡¤Pa/(K¡¤mol)]
temp = 673; %temperature [K]
coninit = zeros(Nx,Ny,Nz)+10000; %Cinit=10000[10^-2ppm];
con = coninit;
dtime=(dx*dx/D)*0.001;

%create matrix for calculating total hydrogen atoms
con_new = zeros(Nx,Ny,Nz); 
dcmatrix  = zeros(Nx,Ny,Nz);
dsmatrix = zeros(Nx,Ny,Nz);
dcsmatrix = zeros(Nx,Ny,Nz);

%FDM method solve differential equations
for istep =1:10000
    for i=1:Nx
	for j=1:Ny
	for k=1:Nz
		jp=j+1;
		jm=j-1;
		ip=i+1;
		im=i-1;
		kp=k+1;
		km=k-1;

		if(im==0)im=Nx;end
		if(ip==(Nx+1))ip=1;end
		if(jm==0)jm = Ny;end
		if(jp==(Ny+1))jp=1;end
		if(km==0)km=Nz;end
		if(kp==(Nz+1))kp=1;end
		
        
		hne=con(ip,j,k); %at (i+1.j,k) "eastern point"
		hnw=con(im,j,k); %at (i-1,j,k) "western point"
		hns=con(i,jm,k); %at (i,j-1,k) "southern point"
		hnn=con(i,jp,k); %at (i,j+1,k) "northern point"
		hnf=con(i,j,kp); %at (i,j,k+1) "upper point"
		hnb=con(i,j,km); %at (i,j,k-1) "under point"
		hnc=con(i,j,k); %at (i,j,k) "centeral point"

		hse=hs(ip,j,k);
		hsw=hs(im,j,k);
		hss=hs(i,jm,k);
		hsn=hs(i,jp,k);
		hsf=hs(i,j,kp);
		hsb=hs(i,j,km);
		hsc=hs(i,j,k);

        dc2dx2=D*(hnw+hne-2.0*hnc)/dx/dx; %dx=dy=dz
        dc2dy2=D*(hnn+hns-2.0*hnc)/dy/dy;
        dc2dz2=D*(hnf+hnb-2.0*hnc)/dz/dz;

        db2dx2=D*(hsw+hse-2.0*hsc)/dx/dx;
        db2dy2=D*(hsn+hss-2.0*hsc)/dy/dy;
        db2dz2=D*(hsf+hsb-2.0*hsc)/dz/dz; 

        de2dx2=D*((hne-hnw)*(hse-hsw))/(4.0*dx*dx);
        de2dy2=D*((hnn-hns)*(hsn-hss))/(4.0*dy*dy);
        de2dz2=D*((hnf-hnb)*(hsf-hsb))/(4.0*dz*dz);

        dcdt=(dc2dx2+dc2dy2+dc2dz2) - v_H/(R*temp)*((de2dx2 + de2dy2 + de2dz2)+con(i,j,k)*(db2dx2+db2dy2+db2dz2));
		dcmatrix(i,j,k) = dc2dx2+dc2dy2+dc2dz2;
        dsmatrix(i,j,k) = v_H/(R*temp)*con(i,j,k)*(db2dx2+db2dy2+db2dz2);
        dcsmatrix(i,j,k) = v_H/(R*temp)*(de2dx2 + de2dy2 + de2dz2);
        
        
        if (con(i,j,k)+dtime*dcdt>=0) 
            con_new(i,j,k)=con(i,j,k)+dtime*dcdt;
        else 
           con_new(i,j,k)=0.0;
        end
           
    end %k
	end %j
	end %i

    for i=1:Nx
    for j=1:Ny
    for k=1:Nz
        con(i,j,k)=con_new(i,j,k);
    end %k
    end %j
    end %i
    
    %output vtk result
	if((mod(istep,nprint)==0) ||(istep==1))
		fprintf('done step: %5d\n',istep);
        vtkmulti('Ctotal.vtk', 'hydrostress', hs,  'InitH', coninit, 'dcmatrix', dcmatrix, 'dsmatrix', dsmatrix, 'dcsmatrix', dcsmatrix, 'Hcon', con)
	end 
end 



function vtkmulti(filename,varargin)
fid = fopen(filename, 'w');
% VTK files contain five major parts
% 1. VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 2.0\n');
% 2. Title
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'ASCII\n');
sx = 1;
sy = 1;
sz = 1;
ox = 0;
oy = 0;
oz = 0;
ndata = length(varargin)/2;
m1 = varargin{2};
[nx, ny, nz] = size(m1);
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
fprintf(fid, 'DIMENSIONS %d %d %d\n', nx, ny, nz);
fprintf(fid, ['SPACING ', num2str(sx), ' ', num2str(sy), ' ', num2str(sz), '\n']);
fprintf(fid, ['ORIGIN ', num2str(ox), ' ', num2str(oy), ' ', num2str(oz), '\n']);
fprintf(fid, 'POINT_DATA %d\n', nx*ny*nz);
for i=1:ndata
    title = varargin{2*i-1};
    m = varargin{2*i};
    fprintf(fid, ['SCALARS ', title, ' double 1\n']);
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%g\t', m);
    fprintf(fid,'\n');
end 

fclose(fid);