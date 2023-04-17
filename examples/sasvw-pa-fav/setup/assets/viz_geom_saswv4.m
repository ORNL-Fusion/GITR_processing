close all
solps_geom = dlmread('solps_rz.txt')';

% SAVE ALL SAVED PLOTS FROM THIS SCRIPT TO '../../plots'

% order indices
A = dlmread('solps_rz.txt');
A = A(1:end-1,:);
ind = zeros(114/2,1);
i = 0;
for j = [1:114]
    if rem(j,2) ~= 0
        i=i+1;
        ind(i) = j;
    end
end
ind(end+1) = 114;

r = A(ind,1)/1000;
z = A(ind,2)/1000;

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