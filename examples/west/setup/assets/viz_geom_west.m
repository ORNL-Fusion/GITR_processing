close all
solps_geom = dlmread('gitr_rz.txt')';

% SAVE ALL SAVED PLOTS FROM THIS SCRIPT TO '../../plots'

% order indices
A = dlmread('gitr_rz.txt');

r = A(:,1);
z = A(:,2);

figure(1);
axis equal
plot(r,z);
hold on
scatter(r,z);
a = [1:length(r)]'; b = num2str(a); c = cellstr(b);
dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
text(r+dx, z+dy, c);
title({'DIII-D SAS-VW4 Cross-Sectional','Geometry'});
xlabel('r [m]')
ylabel('z [m]')
set(gca,'fontsize',16)