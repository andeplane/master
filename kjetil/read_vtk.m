%
% Reads VTK-file format from FEniCS. 
%
function [points, connectivity, offsets, velocities] = read_vtk(file)
fid = fopen(file,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

fgets(fid);   % <?xml version="1.0"?>
fgets(fid);   % <VTKFile type="UnstructuredGrid"  version="0.1"  >
fgets(fid);   % <UnstructuredGrid>
fgets(fid);   % <Piece  NumberOfPoints="XXXX" NumberOfCells="XXXX">
fgets(fid);   % <Points>

points = fgets(fid);

fgets(fid);   % </Points>
fgets(fid);   % <Cells>

connectivity = fgets(fid);
offsets = fgets(fid);
types = fgets(fid);

fgets(fid);   % </Cells>
fgets(fid);   % <PointData  Vectors="u">

velocities = fgets(fid);

points = regexp(points,'(?<=>)(.*\n?)(?=<)','match');
points = str2num(points{1});

connectivity = regexp(connectivity,'(?<=>)(.*\n?)(?=<)','match');
connectivity = str2num(connectivity{1});

offsets = regexp(offsets,'(?<=>)(.*\n?)(?=<)','match');
offsets = str2num(offsets{1});

velocities = regexp(velocities,'(?<=>)(.*\n?)(?=<)','match');
velocities = str2num(velocities{1});
fclose(fid);

end

