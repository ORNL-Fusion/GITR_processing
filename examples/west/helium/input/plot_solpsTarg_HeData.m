data = table2array(readtable('/Users/Alyssa/Dev/GITR/west/helium/output/solpsTarg.txt'));

R_minus_Rsep = data(:,1)
r = data(:,2)
z = data(:,3)
Te = data(:,4)
Ti = data(:,5)
Btot = data(:,12)
Bangle = data(:,13)

Flux = zeros(36,3)
Flux(:,1) = data(:,6)
Flux(:,2) = data(:,7)
Flux(:,3) = data(:,8)

n = zeros(36,3)
n(:,1) = data(:,9)
n(:,2) = data(:,10)
n(:,3) = data(:,11)

figure(1)
hold on
plot(R_minus_Rsep, Te)
plot(R_minus_Rsep, Ti)
title('Temps')
legend({'Te','Ti'})

figure(2)
hold on
plot(R_minus_Rsep, Btot)
title('B field')

figure(3)
hold on
plot(R_minus_Rsep, Bangle)
title('B angle')

figure(4)
hold on
plot(R_minus_Rsep, Flux)
title('Fluxes')
legend({'0','1','2'})

figure(5)
hold on
plot(R_minus_Rsep, n)
title('Densities')
legend({'0','1','2'})

% xlabel('R - R_{sep}') % x-axis label
% ylabel('B [Tesla]') % y-axis label
% set(gca,'fontsize',16)