close all
clear all
do_plot = 1;
filename = 'input/bField.nc'

r = ncread(filename,'r');
z = ncread(filename,'z');

br = ncread(filename,'br');
bt = ncread(filename,'bt');
bz = ncread(filename,'bz');

figure
h = pcolor(r,z,br')
h.EdgeColor = 'none';
colorbar

figure
h = pcolor(r,z,bt')
h.EdgeColor = 'none';
colorbar

figure
h = pcolor(r,z,bz')
h.EdgeColor = 'none';
colorbar

filename = 'input/plasmaProfiles.nc'

gridr = ncread(filename,'gridr');
gridz = ncread(filename,'gridz');
% 
% 
% te = ncread(filename,'te');
%     ret = pcolor_plot(do_plot,gridr,gridz,te,'r [m]','z [m]',{'T_e [eV]','He Discharge'})
ti = ncread(filename,'ti');
    ret = pcolor_plot(do_plot,gridr,gridz,ti,'r [m]','z [m]',{'T_i [eV]','He Discharge'})
% ne = ncread(filename,'ne');
%     ret = pcolor_plot(do_plot,gridr,gridz,ne,'r [m]','z [m]',{'n_e [m^{-3}]','He Discharge'})
ni = ncread(filename,'ni');
    ret = pcolor_plot(do_plot,gridr,gridz,ni,'r [m]','z [m]',{'n_i [m^{-3}]','He Discharge'})
% flux = ncread(filename,'flux');
% mass = ncread(filename,'mass');
%     ret = pcolor_plot(do_plot,gridr,gridz,mass,'r [m]','z [m]',{'Mass [amu]','He Discharge'})
% charge = ncread(filename,'charge');
%     ret = pcolor_plot(do_plot,gridr,gridz,charge,'r [m]','z [m]',{'Charge [#]','He Discharge'})
Br = ncread(filename,'Br');
    ret = pcolor_plot(do_plot,gridr,gridz,Br,'r [m]','z [m]',{'B_r [T]','He Discharge'})
Bt = ncread(filename,'Bt');
    ret = pcolor_plot(do_plot,gridr,gridz,Bt,'r [m]','z [m]',{'B_t [T]','He Discharge'})
Bz = ncread(filename,'Bz');
    ret = pcolor_plot(do_plot,gridr,gridz,Bz,'r [m]','z [m]',{'B_z [T]','He Discharge'})
% 
% Bmag = ncread(filename,'Bmag');
vr = ncread(filename,'vr');
    ret = pcolor_plot(do_plot,gridr,gridz,vr,'r [m]','z [m]',{'v_r [m/s]','He Discharge'})
vt = ncread(filename,'vt');
    ret = pcolor_plot(do_plot,gridr,gridz,vt,'r [m]','z [m]',{'v_t [m/s]','He Discharge'})
vz = ncread(filename,'vz');
    ret = pcolor_plot(do_plot,gridr,gridz,vz,'r [m]','z [m]',{'v_z [m/s]','He Discharge'})
Er = ncread(filename,'Er');
Er(find(Er >0)) = nan;
    ret = pcolor_plot(do_plot,gridr,gridz,Er,'r [m]','z [m]',{'E_r [V/m]','He Discharge'})
Ez = ncread(filename,'Ez');
Ez(find(Ez >0)) = nan;
    ret = pcolor_plot(do_plot,gridr,gridz,Ez,'r [m]','z [m]',{'E_z [V/m]','He Discharge'})

% gradTer = ncread(filename,'gradTer');
%     ret = pcolor_plot(do_plot,gridr,gridz,gradTer,'r [m]','z [m]',{'B dot \nabla T_{er} [eV/m]','He Discharge'})
%     hl = title({'$(\hat{B} \cdot{} \nabla T_{e}) \hat{B}_r [eV/m]$','He Discharge'})
%     set(hl, 'Interpreter', 'latex');
% gradTe = ncread(filename,'gradTe');
%     ret = pcolor_plot(do_plot,gridr,gridz,gradTe,'r [m]','z [m]',{'\hat{B} \dot{} \nabla T_{e} [eV/m]','He Discharge'})
%     hl = title({'$\hat{B} \cdot{} \nabla T_{e} [eV/m]$','He Discharge'})
%     set(hl, 'Interpreter', 'latex');
% gradTet = ncread(filename,'gradTet');
%     ret = pcolor_plot(do_plot,gridr,gridz,gradTet,'r [m]','z [m]',{'B dot \nabla T_{et} [eV/m]','He Discharge'})
%     hl = title({'$(\hat{B} \cdot{} \nabla T_{e}) \hat{B}_t [eV/m]$','He Discharge'})
%     set(hl, 'Interpreter', 'latex');
% gradTez = ncread(filename,'gradTez');
%     ret = pcolor_plot(do_plot,gridr,gridz,gradTez,'r [m]','z [m]',{'B dot \nabla T_{ez} [eV/m]','He Discharge'})
%         hl = title({'$(\hat{B} \cdot{} \nabla T_{e}) \hat{B}_z [eV/m]$','He Discharge'})
%     set(hl, 'Interpreter', 'latex');
gradTi = ncread(filename,'gradTi');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTi,'r [m]','z [m]',{'\nabla T_{i} [eV/m]','He Discharge'})
    hl = title({'$\hat{B} \cdot{} \nabla T_{i} [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTir = ncread(filename,'gradTir');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTir,'r [m]','z [m]',{'B dot \nabla T_{ir} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{i}) \hat{B}_r [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTit = ncread(filename,'gradTit');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTit,'r [m]','z [m]',{'B dot \nabla T_{it} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{i}) \hat{B}_t [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');
gradTiz = ncread(filename,'gradTiz');
    ret = pcolor_plot(do_plot,gridr,gridz,gradTiz,'r [m]','z [m]',{'B dot \nabla T_{iz} [eV/m]','He Discharge'})
        hl = title({'$(\hat{B} \cdot{} \nabla T_{i}) \hat{B}_z [eV/m]$','He Discharge'})
    set(hl, 'Interpreter', 'latex');

% figure(1)
% h = pcolor(gridr,gridz,te')
% h.EdgeColor = 'none';
% colorbar
% set(gca, 'ColorScale', 'log')

Br = -Br;
Bz = -Bz;
Bt = -Bt;
Bnorm = sqrt(Br.^2 + Bt.^2 + Bz.^2);
var = Bt'./Bnorm'.*gradTi';
var_pos = var;
var_neg = var;
var_pos(find(var_pos > 0)) = nan;
var_neg(find(var_neg < 0)) = nan;

figure
 h = pcolor(gridr,gridz,var_pos)
 h.EdgeColor = 'none';
colorbar
% set(gca, 'ColorScale', 'log')

figure
 h = pcolor(gridr,gridz,var_neg)
 h.EdgeColor = 'none';
colorbar
if (exist('x1') == 0)
    fid = fopen('input/gitrGeometry.cfg');

    tline = fgetl(fid);
    tline = fgetl(fid);
    for i=1:10
        tline = fgetl(fid);
        evalc(tline);
    end
    length_side = length;
    clear length;
    Zsurface = Z;
end
%% calculate forces
nu_s = fluid_momentum_loss_rate(ti,184,2,ni,1,1);
    figure
h = pcolor(gridr,gridz,nu_s')
h.EdgeColor = 'none';
colorbar
set(gca,'ColorScale','log')

drag_force = 184*1.67e-27*vt.*nu_s;
h = pcolor(gridr,gridz,drag_force')
h.EdgeColor = 'none';
colorbar
B = beta_i(184,1,2);
thermal_force = B*1.602e-19*gradTit;
figure
h = pcolor(gridr,gridz,thermal_force')
h.EdgeColor = 'none';
colorbar
hold on
for i=1:length(x1)
    plot([x1(i) x2(i)],[z1(i) z2(i)],'k')
%     text(x1(i)+0.1,z1(i)+0.1,string(i))
end
function ret = line_plot(do_plot,x,y,xl,yl,t)
if do_plot
    figure
    plot(x,y,'LineWidth',2)
    title(t)
    ylabel(yl)
    xlabel(xl)
    set(gca,'FontSize',14)
    ret = 1;
else
    ret = 0;
end
end

function ret = pcolor_plot(do_plot,r,z,y,xl,yl,t)
if do_plot
    figure
h = pcolor(r,z,y')
h.EdgeColor = 'none';
colorbar

    title(t)
    ylabel(yl)
    xlabel(xl)
    set(gca,'FontSize',14)
%     axis equal
    axis([min(r) max(r) min(z) max(z)])
    ret = 1;
else
    ret = 0;
end
end
function nu_0 = fluid_momentum_loss_rate(Ti,m,mbackground,n,z,zbackground)
ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;
mbackground = mbackground*MI;
m = m*MI;
vTh = sqrt(2*Ti*Q/mbackground + 2*Ti*Q/m);

lam_d = sqrt(EPS0*Ti./(n*Q));%only one q in order to convert to J
lam = 12*pi*n.*lam_d.^3;
fac = sqrt(m*(Ti*Q).^3);
fac = (1+m/mbackground)^(3/2)*sqrt(m)*(Ti*Q).^(3/2);
% nu_00 = 4*sqrt(2*pi)* Q^4*z^2*(zbackground.^2)*log(lam)*n/(fac*4*4*3*pi*pi*EPS0*EPS0);
% nu_friedberg = sqrt(2)/12/pi^(3/2)* Q^4*z^2*(zbackground.^2)*log(lam)*n/(fac*EPS0*EPS0);
% 
reduced_mass = m*mbackground/(m + mbackground);
% vTh = sqrt(2*Ti*Q/mbackground );
% nu_00 = 4/3/sqrt(pi) * m/reduced_mass * 4*pi*n*Q^4*z^2*zbackground.^2*log(lam)/((4*pi*EPS0)^2*m^2*vTh^3);
vTh = sqrt(2*Ti.*Q./mbackground + 2*Ti.*Q./m);
% vTh = sqrt( 2*Ti*Q/m);
c1 = 0.85*4/3/sqrt(pi) /reduced_mass * 4*pi*Q^4.*z^2.*zbackground.^2/((4*pi*EPS0)^2*m);
nu_0 = c1*log(lam).*n./(vTh.^3);
% nu_0 = 0.85*4/3/sqrt(pi) /reduced_mass * 4*pi*Q^4.*z^2.*zbackground.^2.*log(lam).*n./((4*pi*EPS0)^2*m*vTh.^3);

% nu_0 = 7.35898e5/(2*2*5*1000);
% nu_0 = 8.25e5/2e4;
end
function B = beta_i(mi,Z,mD)
mu = mi/(mD+mi);
B = 3*(mu+5*sqrt(2)*Z^2*(1.1*mu^(5/2) - 0.35*mu^(3/2)) - 1)/(2.6-2*mu + 5.4*mu^2);

end