%function writeMeshToFile(MESH)
fclose all;
NODES = gcoord;
ELEMS = Elem2Node;


tic
filename = 'meshtofile.xml';
fid=fopen(filename,'W');

fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?> \n \n');
fprintf(fid,'<dolfin xmlns:dolfin="http://www.fenics.org/dolfin/"> \n');
fprintf(fid,['<mesh celltype="triangle" dim="' num2str(size(NODES,1)) '"> \n']);
fprintf(fid,['<vertices size="' num2str(size(NODES,2)) '"> \n']);
fprintf(fid,'<vertex index="%d" x="%.16f" y="%.16f"/> \n',[0:(size(NODES,2)-1); NODES(1,:);NODES(2,:)]);
fprintf(fid,'</vertices> \n');
fprintf(fid,['<cells size="' num2str(size(ELEMS,2)) '"> \n']);
fprintf(fid,'<triangle index="%d" v0="%d" v1="%d" v2 ="%d"/> \n',[0:(size(ELEMS,2)-1)';ELEMS(1,:)-1;ELEMS(2,:)-1;ELEMS(3,:)-1]);
fprintf(fid,'</cells> \n');


fprintf(fid,'</mesh> \n');

fprintf(fid,'<mesh_function> \n');
fprintf(fid,['<mesh_value_collection name="boundaries" type="uint" dim="1" size="' num2str(size(FacetID,1)) '"> \n']);
fprintf(fid,'<value cell_index="%d" local_entity="%d" value="%d" /> \n',FacetID');
fprintf(fid,'</mesh_value_collection> \n');
fprintf(fid,'</mesh_function> \n');

fprintf(fid,'</dolfin> \n');




status = fclose(fid);
display(['Mesh written to ' filename ' in ' num2str(toc) 's. File closed with exit status ' num2str(status) '.'])
%end