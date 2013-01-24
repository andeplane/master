function [gcoord,Elem2Node,bc_node,phase_id] = triangle_mesh_generation(coord,segm,hole,regions,nodes,extraoptions,maxArea,foldername)


modelname  =  [foldername 'domain'];
nb_coord   = size(coord,2);
nb_segm    = size(segm,2);
nb_holes   = size(hole,2);
nb_regions = size(regions,2);

%%% Writing .poly file
fid=fopen([modelname,'.poly'],'w');
fprintf(fid,'# %s \n',modelname);
fprintf(fid,'%g %g %g %g \n',nb_coord,2,0,0);
fprintf(fid,'# nodes with coordinates x y \n');
for i=1:nb_coord
    fprintf(fid,'%g %g %g \n',i,coord(:,i));
end
fprintf(fid,'# segments \n');
fprintf(fid,'%g %g \n',nb_segm,1);
for i=1:nb_segm
    fprintf(fid,'%g %g %g %g\n',i,segm(:,i));
end
fprintf(fid,'# holes \n');
fprintf(fid,'%g \n',nb_holes);
for i=1:nb_holes
    fprintf(fid,'%g %g %g %g\n',i,hole(:,i));
end
fprintf(fid,'# regions \n');
fprintf(fid,'%g \n',nb_regions);
for i=1:nb_regions
    fprintf(fid,'%g %g %g %g\n',i,regions(:,i));
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%


%%%%  executing triangle

if nodes>=6
    additionalOptions='-o2';
    nnodel=6;
else
    additionalOptions='';
    nnodel=nodes;
end
[status,result]=system(['./triangleMac -I -p -q32 -Aa -a',num2str(maxArea,'%16f'),' ',additionalOptions,' ' ,extraoptions ' ',modelname,'.poly']);
%!//kvant/felles/kurs/fys-geo4510/2011/Mesh/triangle.exe -pAq32a0.03 //kvant/felles/kurs/fys-geo4510/2011//domain.poly
% also works but no parameters possible:
% !triangle -p -a.001 square.poly     
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%  Reading .node .ele
% reading nodes coordinates x y
fid = fopen([modelname,'.node'], 'r','l');
tmp = fscanf(fid,'%f');
nnod=tmp(1);
tmp(1:4)=[];
tmp=reshape(tmp,4,nnod);
bc_node=find(tmp(4,:)~=0);
bc_node=[bc_node;tmp(4,bc_node)];
gcoord=tmp(2:3,:);
fclose(fid);

% reading elements nodes
fid = fopen([modelname,'.ele'], 'r','l');
tmp = fscanf(fid,'%f');
nel=tmp(1);
tmp(1:3)=[];
tmp=reshape(tmp,2+nnodel,nel); 
Elem2Node=tmp(2:1+nnodel,:)';
phase_id = tmp(end,:);
fclose(fid);
Elem2Node=Elem2Node';
if nodes==7
   nnod=max(max(Elem2Node)); 
   nel=size(Elem2Node,2);
   Elem2Node=[Elem2Node;((1:nel)+nnod)];
   gcx=mean(reshape(gcoord(1,Elem2Node(1:3,:)),3,nel));
   gcy=mean(reshape(gcoord(2,Elem2Node(1:3,:)),3,nel));
   gc=[gcx;gcy];
   gcoord=[gcoord gc];
end

% % plot
% trimesh(Elem2Node,gcoord(1,:),gcoord(2,:),0*gcoord(1,:))
% view(2)
% axis equal
% axis off
% caxis([0 1])
end