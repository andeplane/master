function [dt, points, fn] = triangulateObjFile(filename)

fid = fopen(filename,'r');

tline = fgetl(fid);
points = zeros(1,3);
numPoints = 0;
while ischar(tline) 
    if(numel(tline) > 2 && strcmp(tline(1:3),'v  '))
        % Only use the numbers that are the vector representation starting with
        % 'v  %f %f %f'
        numPoints = numPoints+1;
        points(numPoints,:) = str2num(tline(3:length(tline)));
    end
    
    tline = fgetl(fid);
end
fclose(fid);
points = points*1e-9;

dt = DelaunayTri(points);

[tri Xb] = freeBoundary(dt);
tr = TriRep(tri, Xb);
P = incenters(tr);
fn = faceNormals(tr);
end

