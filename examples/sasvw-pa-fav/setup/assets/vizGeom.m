close all
% SAVE ALL SAVED PLOTS FROM THIS SCRIPT TO '../../plots'

% order indices
A = dlmread('gitr_rz.txt');
A = A(1:end-1,:);
stop = length(A);
ind = zeros(int8(floor(stop))/2,1);
i = 0;
for j = [1:stop]
    %if rem(j,2) ~= 0
        i=i+1;
        ind(i) = j;
   % end
end
ind(end+1) = stop;

r = A(ind,1);%/1000;
z = A(ind,2);%/1000;

figure(1);
axis equal
plot(r,z);
hold on
scatter(r,z);
a = [1:length(r)]'; b = num2str(a); c = cellstr(b);
dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
text(r+dx, z+dy, c, 'fontsize', 15);
title({'DIII-D SAS-VW5 Cross-Sectional','Geometry'});
xlabel('r [m]')
ylabel('z [m]')
axis equal
set(gca,'fontsize',16)