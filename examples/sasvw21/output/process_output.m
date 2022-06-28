fileID = fopen('../setup/assets/geom-SASV/gitr_rz.txt','r');
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
figure(6)
plot(R,Z)
%plot3(R,T,Z)
hold on

for i=1:1:length(x0)
%plot3(x(:,i),y(:,i),z(:,i))
plot(x(:,i),z(:,i))
end

axis equal
%xlim([1.43 1.56])
%ylim([1.07 1.25])
xlabel('r [m]')
zlabel('z [m]')
title('Histories')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat(pwd,'/surface.nc');
grossDep0 = ncread(file,'grossDeposition');
grossEro0 = ncread(file,'grossErosion');

% subset = find(Zsurface);
% x = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
% y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
% z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];



















