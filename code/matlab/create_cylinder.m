voxels = 100;
N = [voxels voxels voxels];
A = zeros(N);
cylinder_radius = 40;
center = 50;

for i=1:voxels
    for j=1:voxels
        for k=1:voxels
            di = i - center;
            dj = j - center;
            dr2 = di*di + dj*dj;
            if dr2 > cylinder_radius^2
                A(i,j,k) = 1;
            end
        end
    end
end

fid = fopen('/projects/master/code/worlds/cylinder_m.bin', 'w');
fwrite(fid, N, 'unsigned char');
fwrite(fid, A, 'unsigned char');
fclose(fid);