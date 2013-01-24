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
fprintf(fid,'<vertex index="%d" x="%.2f" y="%.2f"/> \n',[0:(size(NODES,2)-1)';NODES(1,:);NODES(2,:)]);

% for i = 1:size(NODES,2)
%     fprintf(fid,['<vertex index="' num2str(i-1) '" x="' num2str(NODES(1,i)) '" y="' num2str(NODES(2,i)) '"/> \n']);
% end
fprintf(fid,'</vertices> \n');
fprintf(fid,['<cells size="' num2str(size(ELEMS,2)) '"> \n']);

fprintf(fid,'<triangle index="%d" v0="%d" v1="%d" v2 ="%d"/> \n',[0:(size(ELEMS,2)-1)';ELEMS(1,:)-1;ELEMS(2,:)-1;ELEMS(3,:)-1]);


%fprintf(fid,'<triangle index="0" v0="0" v1="1" v2="3" v3="4"/> \n');
% ...
fprintf(fid,'</cells> \n');
fprintf(fid,'</mesh> \n');
fprintf(fid,'</dolfin> \n');

% 
% <?xml version="1.0" encoding="UTF-8"?>
% 
% <dolfin xmlns:dolfin="http://www.fenics.org/dolfin/">
%   <mesh celltype="tetrahedron" dim="3">
%     <vertices size="8">
%       <vertex index="0" x="0" y="0" z="0" />
%    <vertex index="7" x="0" y="1" z="1" />
%     </vertices>
%     <cells size="6">
%       <tetrahedron index="0" v0="0" v1="1" v2="3" v3="4"/>
% 
%             <tetrahedron index="5" v0="1" v1="3" v2="6" v3="7"/>
%     </cells>
%   </mesh>
% </dolfin>




status = fclose(fid);
display(['Mesh written to ' filename ' in ' num2str(toc) 's. File closed with exit status ' num2str(status) '.'])
%end