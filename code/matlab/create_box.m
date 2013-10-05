% Creates a box with 
%
%
%
%
function create_box(fraction)
    voxels = 100;
    N = [voxels voxels voxels];
    A = zeros(N);
    center = 50;

    for i=1:voxels
        for j=1:voxels
            for k=1:voxels
                dj = j - center;

                if abs(dj)>center*fraction
                    A(i,j,k) = 1;
                end
            end
        end
    end
    
    fid = fopen('/projects/master/code/worlds/box_fraction_0.2_m.bin', 'w');
    fwrite(fid, N, 'uint');
    fwrite(fid, A, 'unsigned char');
    fclose(fid);
end