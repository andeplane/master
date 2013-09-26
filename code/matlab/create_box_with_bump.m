% Creates a box with a slice
function create_box_with_bump()
    voxels = 100;
    N = [voxels voxels voxels];
    A = zeros(N);
    center = voxels / 2;

    for i=1:voxels
        for j=1:voxels
            for k=1:voxels
                dj = j - center;

                if abs(dj)>center*0.2
                    A(i,j,k) = 1;
                end
                
                if dj>(center*0.2) && abs(dj)<(center*0.2+3) && abs(i-center) < 2
                    A(i,j,k) = 0;
                end
            end
        end
    end
    
    x = zeros(voxels,1);
    y = zeros(voxels,1);
    z = zeros(voxels,1);
    count = 1;
    for i=1:voxels
        for j=1:voxels
            for k=1:voxels
                if A(i,j,k)==1
                    x(count) = i;
                    y(count) = j;
                    z(count) = k;
                    count = count + 1;
                end
            end
        end
    end

    %plot3(x,y,z,'ko')
    plot(x,y,'ko')
    axis('equal')
    
    fid = fopen('/projects/master/code/worlds/box_with_bump_m.bin', 'w');
    fwrite(fid, N, 'uint');
    fwrite(fid, A, 'unsigned char');
    fclose(fid);
end