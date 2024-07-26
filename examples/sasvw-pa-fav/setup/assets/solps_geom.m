%plot org

clear all
A = readmatrix("sas-vw_v005_mod.txt")./1000;

%%
figure(1)
hold on
for i=1:length(A)-1
    ogr = plot(A(i:i+1,1),A(i:i+1,2),'-ob','Linewidth',1,'HandleVisibility','off');
end
%plot([A(end,1),A(1,1)],[A(end,2),A(1,2)],'-om','DisplayName','closing point');
%print 'plotted ogr'

%% plot corners and topcut

nxny =  read_b2f_variable("b2fgmtry", "nx,ny");

crx_long =  read_b2f_variable("b2fgmtry", "crx");
cry_long =  read_b2f_variable("b2fgmtry", "cry");

nx = nxny(1);
ny = nxny(2);
crx0= reshape(crx_long,nx+2,ny+2,4);
cry0= reshape(cry_long,nx+2,ny+2,4);

crx = crx0;
cry = cry0;
crx(:,:,4) = crx0(:,:,3);
cry(:,:,4) = cry0(:,:,3);
crx(:,:,3) = crx0(:,:,4);
cry(:,:,3) = cry0(:,:,4);
%print 'got corners'
%print 'plotting b2fgmtry'

%% plot GITR grid
profilesFile = strcat('../../input/plasmaProfiles.nc');
r = ncread(profilesFile,'gridr');
z = ncread(profilesFile,'gridz');
nr = length(r);
nz = length(z);
dr = r(4)-r(3);
dz = z(4)-z(3);

figure(10)
xticks(1.4:dr:1.52);
yticks(1.08:dz:1.24);
ax = gca;
ax.GridColor = [0 0 0];
grid on
hold on

%%
%figure(2)
hold on
for i=1:nx+2
    for j=1:ny+2
        b2f = patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),nan); %, 'DisplayName','b2fgmtry corners')
    end
end
%print 'plotted b2fgmtry'


%get from gitrGeometry.cfg

fid = fopen('../../input/gitrGeometry.cfg');

tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

W_indices = [70:70+110];
x1 = [];
z1 = [];

for i = 5:241
    x1 = [x1, str2double(strip(fgetl(fid)))];
end

tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

for i = 245:481
    z1 = [z1, str2double(strip(fgetl(fid)))];
end

%set real wall coords for adjusted 3rd leg
r_real = [1502, 1485.5, 1491.05, 1499.909817]/1000;
z_real = [1163, 1123, 1112.75, 1110.921234]/1000;

hold on
gitr = plot(x1,z1,'b','linewidth',1);
W = plot(x1(W_indices),z1(W_indices),'m','linewidth',1);
scatter(x1(W_indices),z1(W_indices),50,'.','m');
%real_leg = plot(r_real, z_real, 'color', "#77AC30", 'linewidth',2);
%strikepoint = plot([1.49829829],[1.19672716],'pentagram');
%strikepoint.MarkerFaceColor = [1 0.5 0];
%strikepoint.MarkerSize = 20;
%legend([ogr,b2f,gitr], 'ogr wall','b2fgmtry target','GITR boundary')
[h, icons] = legend([gitr,W],'Carbon','Tungsten','fontsize',40);
%[h, icons] = legend([gitr,W,real_leg,strikepoint],'Carbon','Tungsten','Real W 3rd Leg','Strikepoint','fontsize',40,'Location','northwest');

axis equal
%title({'Cross Section of SAS-VW Geometry', 'and SOLPS-ITER Grid'},'fontsize',40)
title('Cross Section of SOLPS-ITER Grid','fontsize',40)
xlabel('R [m]','fontsize',40)
ylabel('Z [m]','fontsize',40)
ax = gca;
ax.FontSize = 30; 
xlim([0.9, 2.5]);
ylim([-1.4, 1.46]);
%xlim([1.39, 1.51]);
%ylim([1.09, 1.24]);


%%
topcut =  read_b2f_variable("b2fgmtry", "topcut")+2;
hold on
figure(1)
j = topcut+15;
for i=1:nx+2
        patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),'blue') %,'DisplayName','b2fgmtry corners with topcut')
end
print 'plotted topcut'

%%
figure(3)
hold on
plot(crx(1:5,topcut,1),cry(1:5,topcut,1),'r','LineWidth',4,'DisplayName','strikepoint')
plot(crx(1:5,topcut,2),cry(1:5,topcut,2),'g','LineWidth',2,'DisplayName','strikepoint')
%print 'plotted strikepoint'
axis equal

%% define rmrs_edges
r_edges = [crx(1,:,3) crx(1,end,4)];
z_edges = [cry(1,:,3) cry(1,end,4)];
scatter(r_edges,z_edges,'r')

for i=1:length(r_edges)
    text(r_edges(i),z_edges(i),string(i)) % scatter points along target edge
end

segment_length = zeros(length(r_edges)-1,1);
r_center = 0*segment_length;
z_center = 0*r_center;
for i=1:length(r_edges)-1
    segment_length(i) = sqrt((r_edges(i) - r_edges(i+1))^2 + (z_edges(i) - z_edges(i+1))^2);
    r_center(i) = 0.5*(r_edges(i) + r_edges(i+1));
    z_center(i) = 0.5*(z_edges(i) + z_edges(i+1));
end

rmrs_edges = [0; cumsum(segment_length)] - sum(segment_length(1:topcut-1));

% plot parallel velocity over the b2fgmtry grid

te =  read_b2f_variable("b2fstate", "te");
ti =  read_b2f_variable("b2fstate", "ti");
gradTit = read_b2f_variable("b2fstate", "gradTit");
bb =  read_b2f_variable("b2fgmtry", "bb");
bb = reshape(bb,nx+2,ny+2,4);

ua =  read_b2f_variable("b2fstate", "ua"); % parallel velocity
ua = reshape(ua,nx+2,ny+2,9);

na =  read_b2f_variable("b2fstate", "na"); 
na = reshape(na,nx+2,ny+2,9); 
naD0 =  na(:,:,1);  % D0 density
naD1 = na(:,:,2); % D1+ density
naC0 = na(:,:,3); % C0 density
naC1 = na(:,:,4); % C1+ density
naC2 = na(:,:,5); % C2+ density
naC3 = na(:,:,6); % C3+ density
naC4 = na(:,:,7); % C4+ density
naC5 = na(:,:,8); % C5+ density
naC6 = na(:,:,9); % C6+ density

%carbon fraction
naC = naC0 + naC1 + naC2 + naC3 + naC4 + naC5 + naC6;
ni = naD0 + naD1 + naC
Cfraction = 100 .* naC./ni;

%figure(4)
for i=1:nx+2
    for j=1:ny+2
        patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),Cfraction(i,j))
    end
end

%% plot te over the b2fgmtry grid

gradTi = reshape(gradTi,nx+2,ny+2);
% ti = reshape(ti,nx+2,ny+2);
figure(5)
for i=1:nx+2
    for j=1:ny+2
        patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),te(i,j))
    end
end

hold on
scatter(r_center,z_center) % plot r_mid along the target edge

%%
rmrs_center = rmrs_edges(1:end-1) + 0.5*segment_length;
te_cells1 = te(1,:);
te_cells2 = te(2,:);
hx =  read_b2f_variable("b2fgmtry", "hx");
hx = reshape(hx,nx+2,ny+2);

hx_cells1 = hx(1,:);
hx_cells2 = hx(2,:);
te_surface = (te_cells1.*0.5.*hx_cells2 + te_cells2.*0.5.*hx_cells1)./(0.5*(hx_cells2 + hx_cells1));

figure(6)
plot(rmrs_center,te_surface)
%%
function field =  read_b2f_variable(solps_geometry_filename, field_name)
%     f = fopen(solps_geometry_filename, 'r')
txt = readlines(solps_geometry_filename);
%     f.close()
%
field_start = 0;
field_end = 0;
found = 0;
%
%      for count, line in enumerate(txt):
for i=1:length(txt)
    %         if found == 0:
    if found == 0
        %             if '*cf' in line:
        if strfind(txt(i),"*cf")
            i;
            %                 words = line.split()
            words = strsplit(txt(i));
            %                 if words[-1] == field_name:
            if strcmp(words(end-1),field_name)
                i;
                field_start = i+1
                found = 1;
            end
        end
    elseif found == 1
        %             if '*cf' in line:
        if strfind(txt(i),"*cf")
            field_end = i;
            found = 2;
        end
    elseif found == 2
        break
    end
end

%
field = [];
txt_list = txt(field_start:field_end-1);
%     for sublist in txt_list:
for j=1:length(txt_list)
    field_previous = field;
    field = [field_previous;str2num(txt_list(j))'];

    %         split_sublist = sublist.split()
    %         for element in split_sublist:
    %             field.append(element)
    %
    %     field = np.array(field)
    %     field = field.astype(np.float)

end
field = reshape(field',[],1);

end