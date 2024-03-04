function Css_Cpp()

%load excel crystal plasticity finite element data
A = xlsread('G120x120x10_Ngrains594_random_tensionX_inc300.xlsx','B:D');
B = xlsread('G120x120x10_Ngrains594_random_tensionX_inc300','F:H');
C = [A, B];

%create axis and hydrostatic stress arrays[Pa]
x=C(:,1); %create x-axis array
y=C(:,2); %create y-axis array
z=C(:,3); %create z-axis array
delta0=((C(:,4)+C(:,5)+C(:,6))/3.0)*10^6; %create hydrostatic stress arrays
delta11=C(:,4)*10^6; %create delta11 direction stress arrays
delta22=C(:,5)*10^6; %create delta22 direction stress arrays
delta33=C(:,6)*10^6; %create delta33 direction stress arrays
length(x)

%number of 3D meshes 90*100*10
Nx=max(x)+1;
Ny=max(y)+1;
Nz=max(z)+1;

hs=zeros(Nx,Ny,Nz);
s11=zeros(Nx,Ny,Nz);
s22=zeros(Nx,Ny,Nz);
s33=zeros(Nx,Ny,Nz);

%create a 3D matrix of hydrostatic stress
for i=1:size(C,1) 
    hs(x(i)+1,y(i)+1,z(i)+1)=delta0(i);
    s11(x(i)+1,y(i)+1,z(i)+1)=delta11(i);
    s22(x(i)+1,y(i)+1,z(i)+1)=delta22(i);
    s33(x(i)+1,y(i)+1,z(i)+1)=delta33(i);
end

time0=clock();
format long;

% calculate steps [m]
dx=1.0e-6;
dy=1.0e-6;
dz=1.0e-6;

%material property parameters
nprint=200; %total output
D = 2.88e-10; %hydrogen atom diffusion parameter [m2 s-1]
v_H = 1.67e-6; %patitial molar mass of hydrogen atom [m3 mol-1]
R = 8.314; %gas constant [N m kmol-1K-1] [m3 ¡¤Pa/(K¡¤mol)]
temp = 673; %temperature [K]
TSSd = 100.0; %which depends on temperature
TSSp=200.0; %which depends on temperature
ax=62.3; %[s-1/2]
Qx = 4.12*10e4; %[J/mol]
qrate = 100;

%create matrix for calculating total hydrogen atoms
dcmatrix  = zeros(Nx,Ny,Nz);
dsmatrix = zeros(Nx,Ny,Nz);
dcsmatrix = zeros(Nx,Ny,Nz);
dpmatrix = zeros(Nx,Ny,Nz);
flagmatrix = zeros(Nx,Ny,Nz);

%calculate Css and Cpp seperately
%In the initial state, all hydrogen atoms are dissolved in the zirconium alloy matrix, [con] is Css, Cpp=0
cinit = 10000;
coninit = zeros(Nx,Ny,Nz)+cinit;
css = coninit;
cpp = zeros(Nx,Ny,Nz);
con = css+cpp;
css_new=zeros(Nx,Ny,Nz);
cpp_new=zeros(Nx,Ny,Nz); 
Deff = D;
dtime=(dx*dx/Deff)*0.001;

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

        dc2dx2=Deff*(hnw+hne-2.0*hnc)/dx/dx;
        dc2dy2=Deff*(hnn+hns-2.0*hnc)/dy/dy;
        dc2dz2=Deff*(hnf+hnb-2.0*hnc)/dz/dz; %dx=dy=dz

        db2dx2=Deff*(hsw+hse-2.0*hsc)/dx/dx;
        db2dy2=Deff*(hsn+hss-2.0*hsc)/dy/dy;
        db2dz2=Deff*(hsf+hsb-2.0*hsc)/dz/dz; 

        de2dx2=Deff*((hne-hnw)*(hse-hsw))/(4.0*dx*dx);       
        de2dy2=Deff*((hnn-hns)*(hsn-hss))/(4.0*dy*dy);
        de2dz2=Deff*((hnf-hnb)*(hsf-hsb))/(4.0*dz*dz);

        dcdt=(dc2dx2+dc2dy2+dc2dz2) - v_H/(R*temp)*((de2dx2 + de2dy2 + de2dz2)+con(i,j,k)*(db2dx2+db2dy2+db2dz2));
       
       
        %kinetics of precipitation and dissolution
        if (css(i,j,k) <=TSSd) && (cpp(i,j,k)==0.0) % diffusion only
            css_new(i,j,k) = css(i,j,k) + dtime*dcdt;
            cpp_new(i,j,k) = cpp(i,j,k);
            flagmatrix(i,j,k) = 1;
            dpmatrix(i,j,k) = 0.0;
            Deff = D/(1+dpdt/dcdt);
        elseif (css(i,j,k) <TSSd) && (cpp(i,j,k)>0.0) % dissolve
            kesi = ax * exp(-Qx/R/temp);
            dpdt = -qrate * kesi*kesi*(TSSp - css(i,j,k)); % negative
            if -dpdt >= cpp(i,j,k) dpdt = -cpp(i,j,k); end
            css_new(i,j,k) = css(i,j,k) + (dcdt - dpdt)*dtime;
            cpp_new(i,j,k) = cpp(i,j,k) + dpdt*dtime;
            flagmatrix(i,j,k) = 2;
            dpmatrix(i,j,k) = dpdt;
            Deff = D/(1+dpdt/dcdt);
        elseif (css(i,j,k) >TSSd) && (css(i,j,k) <=TSSp) % hysteresis
            css_new(i,j,k) = css(i,j,k) + dtime*dcdt;
            cpp_new(i,j,k) = cpp(i,j,k);
            flagmatrix(i,j,k) = 3;
            dpmatrix(i,j,k) = 0.0;
            Deff = D/(1+dpdt/dcdt);
        elseif css(i,j,k) >TSSp %precipitation
            kesi = ax * exp(-Qx/R/temp);
            dpdt = kesi*kesi*(css(i,j,k) - TSSp); %positive
            css_new(i,j,k) = css(i,j,k) + (dcdt - dpdt)*dtime;
            cpp_new(i,j,k) = cpp(i,j,k) + dpdt*dtime;
            flagmatrix(i,j,k) = 4;
            dpmatrix(i,j,k) = dpdt;
            Deff = D/(1+dpdt/dcdt);
        end 
        

	end %k
	end %j
	end %i
    css=css_new;
    cpp=cpp_new;
    con=css_new + cpp_new;
    
	%output vtk result
	if((mod(istep,nprint)==0) ||(istep==1))
		
        if istep==1
            a(1,1)=cinit;
        else
            a(istep/200+1,1)=css_new(38,26,1);
        end
        fprintf('done step: %5d\n',istep);
        vtkmulti('Css_Cpp.vtk', 'hydrostress', hs, 'stress11', s11, 'stress22', s22, 'stress33', s33, 'InitH', coninit, 'dcmatrix', dcmatrix, 'dsmatrix', dsmatrix, ...
        'dcsmatrix', dcsmatrix, 'dpmatrix', dpmatrix,'css', css_new, 'cpp', cpp_new, 'flag', flagmatrix);
    end 

end 

compute_time=etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);


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

