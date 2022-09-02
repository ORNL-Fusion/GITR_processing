close all
solps_geom = dlmread('solps_rz.txt')';

% SAVE ALL SAVED PLOTS FROM THIS SCRIPT TO '../../plots'

% order indices
indices = zeros(1,122/2);
i=0
for j = [1:122]
    if rem(j,2) ~= 0
        i=i+1;
        indices(i) = j;
    end
end
indices = [indices,122:160];

r = solps_geom(1,indices);
z = solps_geom(2,indices);
r1 = r(1:end-1);
r2 = r(2:end);
z1 = z(1:end-1);
z2 = z(2:end);

figure(1)
axis equal
plot(r,z)
hold on
scatter(r,z)
a = [1:length(r)]'; b = num2str(a); c = cellstr(b);
dx = 0.01; dy = 0.01; % displacement so the text does not overlay the data points
text(r+dx, z+dy, c);
title({'DIII-D SAS-VW Cross-Sectional','Geometry'})
xlabel('r [m]')
ylabel('z [m]')
set(gca,'fontsize',16)