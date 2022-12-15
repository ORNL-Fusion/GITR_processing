close all
clear all

% Number of simulated particles
nP = 1e6;

radius = 0.0005;

r = radius*sqrt(rand(nP,1));
theta = 2*pi*rand(nP,1);

% Sampling of  x,y,z
x_sample = r .* cos(theta);
y_sample = r .* sin(theta);
z_sample = zeros(nP,1) + 1e-6;

m = 184;

nPoints = 200;
maxE = 20;
Eb = 8.79; % Surface binding energy
a = 5; % Free parameter
E = linspace(0,maxE,nPoints);
dE = E(2);
thompson2 = a*(a-1)*E.*Eb^(a-1)./(E+Eb).^(a+1);
figure(234)
plot(E,thompson2)
xlabel('E [eV]')
ylabel('pdf')

title('Energy Distribution')


ecdf = cumsum(thompson2);
ecdf = ecdf./ecdf(end);
rand1 = rand(1,nP);

randTheta = 2*pi*rand(1,nP);
randPhi = 0.5*pi*rand(1,nP);
Esamp = interp1(ecdf,E,rand1);

v = sqrt(2*Esamp*1.602e-19/m/1.66e-27)';

    vx = v'.*cos(randTheta).*sin(randPhi);
    vy = v'.*sin(randTheta).*sin(randPhi);
    vz = v'.*cos(randPhi);





ncid = netcdf.create(strcat('./particle_source_0p5mm.nc'),'NC_WRITE')
 
dimP = netcdf.defDim(ncid,'nP',nP);

xVar = netcdf.defVar(ncid,'x','double',[dimP]);
yVar = netcdf.defVar(ncid,'y','double',[dimP]);
zVar = netcdf.defVar(ncid,'z','double',[dimP]);
vxVar = netcdf.defVar(ncid,'vx','double',[dimP]);
vyVar = netcdf.defVar(ncid,'vy','double',[dimP]);
vzVar = netcdf.defVar(ncid,'vz','double',[dimP]);

netcdf.endDef(ncid);
 
netcdf.putVar(ncid, xVar, x_sample);
netcdf.putVar(ncid, yVar, y_sample);
netcdf.putVar(ncid, zVar, z_sample);
netcdf.putVar(ncid, vxVar, vx);
netcdf.putVar(ncid, vyVar, vy);
netcdf.putVar(ncid, vzVar, vz);

netcdf.close(ncid);




