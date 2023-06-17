close all
clear all

A = readmatrix("gitr_rz.txt",'NumHeaderLines',1)./1000;
figure
hold on
for i=1:length(A)-1
    plot(A(i:i+1,1),A(i:i+1,2),'-ok');
end
plot([A(end,1),A(1,1)],[A(end,2),A(1,2)],'-ok');

crx_long =  read_b2f_variable("b2fgmtry", "crx");
cry_long =  read_b2f_variable("b2fgmtry", "cry");

nx = 95;
ny = 40;
crx0= reshape(crx_long,nx+2,ny+2,4);
cry0= reshape(cry_long,nx+2,ny+2,4);

crx = crx0;
cry = cry0;
crx(:,:,4) = crx0(:,:,3);
cry(:,:,4) = cry0(:,:,3);
crx(:,:,3) = crx0(:,:,4);
cry(:,:,3) = cry0(:,:,4);
figure(1)
hold on
for i=1:nx+2
    for j=1:ny+2
        patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),0)
    end
end
topcut =  read_b2f_variable("b2fgmtry", "topcut")+2;

figure(1)
hold on
plot(crx(1:5,topcut,1),cry(1:5,topcut,1),'r','LineWidth',2)
plot(crx(1:5,topcut,2),cry(1:5,topcut,2),'g','LineWidth',2)

r_edges = [crx(1,:,3) crx(1,end,4)];
z_edges = [cry(1,:,3) cry(1,end,4)];
scatter(r_edges,z_edges,'r')

for i=1:length(r_edges)
    text(r_edges(i),z_edges(i),string(i))
end


segment_length = zeros(length(r_edges)-1,1);
r_center = 0*segment_length;
z_center = 0*r_center;
for i=1:length(r_edges)-1
    segment_length(i) = sqrt((r_edges(i) - r_edges(i+1))^2 + (z_edges(i) - z_edges(i+1))^2);
    r_center(i) = 0.5*(r_edges(i) + r_edges(i+1));
    z_center(i) = 0.5*(z_edges(i) + z_edges(i+1));
end

rmrs_edges = [ 0; cumsum(segment_length)] - sum(segment_length(1:topcut-1));

te =  read_b2f_variable("b2fstate", "te");
te = reshape(te,nx+2,ny+2);
figure
for i=1:nx+2
    for j=1:ny+2
        patch(reshape(crx(i,j,:),1,4),reshape(cry(i,j,:),1,4),te(i,j))
    end
end

hold on
scatter(r_center,z_center)

rmrs_center = rmrs_edges(1:end-1) + 0.5*segment_length;
te_cells1 = te(1,:);
te_cells2 = te(2,:);
hx =  read_b2f_variable("b2fgmtry", "hx");
hx = reshape(hx,nx+2,ny+2);

hx_cells1 = hx(1,:);
hx_cells2 = hx(2,:);
te_surface = (te_cells1.*0.5.*hx_cells2 + te_cells2.*0.5.*hx_cells1)./(0.5*(hx_cells2 + hx_cells1));
figure
plot(rmrs_center,te_surface)
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
                field_start = i+1;
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
    field = [field;str2num(txt_list(j))];
%     split_sublist = sublist.split();
%     for element = 1:length(split_sublist)
%         field.append(split_sublist(element))
%   
%     field = np.array(field)
%     field = field.astype(np.float)
%     end

end
field = reshape(field',[],1);

end