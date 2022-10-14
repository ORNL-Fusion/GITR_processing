fileID = fopen('../setup/assets/gitr_rz.txt','r');
formatSpec='%f %f';
sizeA = [2 Inf];
geom = fscanf(fileID,formatSpec,sizeA); 
R = geom(1,:);
T = zeros(length(R));
Z = geom(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% startPosition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat(pwd,'/../input/particleSource.nc');
x0 = ncread(file,'x');
y0 = ncread(file,'y');
z0 = ncread(file,'z');

figure(1)
plot(R,Z)
hold on
scatter(x0,z0)

axis equal
xlim([1.46 1.52])
ylim([1.16 1.25])
xlabel('r [m]')
ylabel('z [m]')
title('Start Positions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% endPosition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_tracks = 1;
file = strcat(pwd,'/positions.nc');
hitWall = ncread(file,'hitWall');
nHit = length(find(hitWall));
hasHit = find(hitWall);
notHit = find(hitWall==0);
x0 = ncread(file,'x');
y0 = ncread(file,'y');
z0 = ncread(file,'z');
vx0 = ncread(file,'vx');
vy0 = ncread(file,'vy');
vz0 = ncread(file,'vz');
time0=ncread(file,'time');
distTraveled = ncread(file,'distTraveled');
charge0 = ncread(file,'charge');
weight0 = ncread(file,'weight');
vtot = sqrt(vx0.^2 +vy0.^2 + vz0.^2);
%E = 0.5*184*1.66e-27*vtot.^2/1.602e-19; %0.5 * amu * vtot.^2 * constants
%figure(2)
%histogram(E)

charge = ncread(file,'charge');
charge_avg = mean(charge)

figure(4)
plot(R,Z)
hold on
scatter(x0(hasHit), z0(hasHit), 'g')
hold on
scatter(x0(notHit), z0(notHit), 'r')

axis equal
xlim([1.46 1.52])
ylim([1.16 1.25])
xlabel('r [m]')
ylabel('z [m]')
title('End Positions')
legend('wall','hasHit','notHit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_tracks
file = strcat(pwd,'/history.nc');
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
sizeArray = size(x);
nP = sizeArray(2);
hit = find(weight(end,:) < 1);

%r = sqrt(x.^2 + y.^2);
figure(6) %3D plotting
plot3(R,T,Z, 'LineWidth',2)
hold on

for i=1:1:length(x0)
plot3(x(:,i),y(:,i),z(:,i), 'LineWidth',0.1)
end

axis equal
xlim([1.43 1.56])
ylim([-0.05 0.05])
zlim([1.07 1.25])
xlabel('r [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Histories')

figure(7) %2D plotting
plot(R,Z,'LineWidth',2)
hold on

for i=1:1:length(x0)
plot(x(:,i),z(:,i),'LineWidth',1)
end

axis equal
xlim([1.43 1.56])
ylim([1.07 1.25])
xlabel('r [m]')
ylabel('z [m]')
title('Histories')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat(pwd,'/surface.nc');
grossDep = transpose(ncread(file,'grossDeposition'));
grossEro = transpose(ncread(file,'grossErosion'));
netEro=grossEro-grossDep

%find where the hell the surface is
Rsurf = R(31:45)
Zsurf = Z(31:45)
figure(10)
plot(Rsurf,Zsurf)
axis equal

figure(11)
plot(Zsurf,grossEro,'-r')
hold on
plot(Zsurf,grossDep,'-g')
hold on
plot(Zsurf,netEro,'-b')

xlabel('z [m]')
ylabel('Particles per Second')
legend('Gross Erosion', 'Redposition', 'Net Erosion')
title('GITR Predicted Erosion and Redeposition Profiles')




















