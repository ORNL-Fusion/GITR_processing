close all
d = dlmread('assets/vvfile')';

% order indices
indices = [1:5,7,6,8:110,112:117];
r = d(1,indices);
z = d(2,indices);
r1 = r(1:end-1);
r2 = r(2:end);
z1 = z(1:end-1);
z2 = z(2:end);

figure(1)
title({'DIII-D SAS-VW Cross-Sectional','Geometry'})
xlabel('r [m]')
ylabel('z [m]')
set(gca,'fontsize',16)
axis equal
plot(r,z)
hold on
scatter(r,z)
a = [1:length(r)]'; b = num2str(a); c = cellstr(b);
dx = 0.01; dy = 0.01; % displacement so the text does not overlay the data points
text(r+dx, z+dy, c);