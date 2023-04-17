fileID = fopen('../../setup/assets/gitr_rz.txt','r');
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
xlim([min(R) max(R)])
ylim([min(Z) max(Z)])
xlabel('r [m]')
ylabel('z [m]')
title('Start Positions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% endPosition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
xlim([min(R) max(R)])
ylim([min(Z) max(Z)])
xlabel('r [m]')
ylabel('z [m]')
title('End Positions')
legend('wall','hasHit','notHit')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_tracks = 1;
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

    figure(6) %3D plotting
    plot3(R,T,Z, 'LineWidth',2)
    hold on

    for i=1:1:length(x0)
    plot3(x(:,i),y(:,i),z(:,i), 'LineWidth',0.1)
    end

    axis equal
    xlim([min(R) max(R)])
    ylim([-0.01 0.01])
    zlim([min(Z) max(Z)])
    xlabel('r [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    title('Histories')

    figure(7) %2D plotting
    plot(R,Z,'LineWidth',2)
    hold on
    
    for i=1:1:length(x(0))
        if charge(i)==1
            plot(x(:,i),z(:,i),'LineWidth',1,'color','red')
        else if charge(i)==2
            plot(x(:,i),z(:,i),'LineWidth',1,'color','blue')
        else if charge(i)==3
            plot(x(:,i),z(:,i),'LineWidth',1,'color','yellow')
        else if charge(i)==4
            plot(x(:,i),z(:,i),'LineWidth',1,'color','green')
        else if charge(i)==5
            plot(x(:,i),z(:,i),'LineWidth',1,'color','cyan')
        else if charge(i)==0
            plot(x(:,i),z(:,i),'LineWidth',1,'color','black')
        else
            plot(x(:,i),z(:,i),'LineWidth',1,'color','magenta')
            end
            end
            end
            end
            end
        end
    end

    axis equal
    xlim([min(R) max(R)])
    ylim([min(Z) max(Z)])
    xlabel('r [m]')
    ylabel('z [m]')
    title('Histories')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% history segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8) %2D plotting
plot(R,Z,'LineWidth',2)
hold on

    
    for i=1:1:length(x)
        if z(1,i) >= 1.125 && z(1,i) <= 1.13
            if charge(i)==1
                plot(x(:,i),z(:,i),'LineWidth',1,'color','red')
            else if charge(i)==2
                plot(x(:,i),z(:,i),'LineWidth',1,'color','blue')
            else if charge(i)==3
                plot(x(:,i),z(:,i),'LineWidth',1,'color','yellow')
            else if charge(i)==4
                plot(x(:,i),z(:,i),'LineWidth',1,'color','green')
            else if charge(i)==5
                plot(x(:,i),z(:,i),'LineWidth',1,'color','cyan')
            else if charge(i)==0
                plot(x(:,i),z(:,i),'LineWidth',1,'color','black')
            else
                plot(x(:,i),z(:,i),'LineWidth',1,'color','magenta')
                end
                end
                end
                end
                end
            end
        end
    end
        axis equal
        xlim([min(R) max(R)])
        ylim([min(Z) max(Z)])
        xlabel('r [m]')
        ylabel('z [m]')
        title('Histories')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = strcat(pwd,'/surface_nP1e5_nT1e6.nc');
grossDep = flip(ncread(file,'grossDeposition'));
grossEro = flip(ncread(file,'grossErosion'));
netEro=grossEro-grossDep;

%find where the hell the surface is
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
legend('Gross Erosion', 'Redposition', 'Net Erosion','Location','northwest')
title('GITR Predicted Erosion and Redeposition Profiles')




















