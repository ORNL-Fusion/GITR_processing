%plot org

clear all
A = readmatrix("sas-vw_v005_mod.txt")./1000;

%%
figure(1)
hold on
for i=1:length(A)-1
    ogr = plot(A(i:i+1,1),A(i:i+1,2),'-oc','Linewidth',1,'HandleVisibility','off');
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
%%
%figure(2)
hold on
for i=1:nx+2
    for j=1:ny+2
        b2f = patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),nan); %, 'DisplayName','b2fgmtry corners')
    end
end
%print 'plotted b2fgmtry'

%%
%get from gitrGeometry.cfg

fid = fopen('gitrGeometry.cfg');

tline = fgetl(fid);

tline = fgetl(fid);

for i=1:2

    tline = fgetl(fid);

    evalc(tline);

end
hold on
gitr = plot(x1,z1,'m');

legend([ogr,b2f,gitr], 'ogr wall','b2fgmtry target','GITR boundary')
axis equal

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
bb =  read_b2f_variable("b2fgmtry", "bb");
bb = reshape(bb,nx+2,ny+2,4);

ua =  read_b2f_variable("b2fstate", "ua"); % parallel velocity
ua = reshape(ua,nx+2,ny+2,9);
%figure(4)
for i=1:nx+2
    for j=1:ny+2
        patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),ua(i,j,2))
    end
end

%% plot te over the b2fgmtry grid

te = reshape(te,nx+2,ny+2);
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