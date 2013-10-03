% Creates a box with a slice
function create_box_with_rotation_2(theta, voxels)
    N = [voxels voxels voxels];
    A = zeros(N);
    center = round(voxels / 2);
    a = [center center];
    n = [sin(theta) cos(theta)];
    
    for i=1:voxels
        i
        for j=1:voxels
            for k=1:voxels
                p = [i j];
                
                % Check nearest distance to line in the channel
                best_dist = norm( (a-p) - ((a-p)*n')*n);
                
                % Check periodic boundaries in x-direction
                p(1) = p(1) + voxels;
                dist = norm( (a-p) - ((a-p)*n')*n);
                if(dist < best_dist)
                    best_dist = dist;
                end
                
                % Check periodic boundaries in y-direction
                p(1) = p(1) - voxels;
                p(2) = p(2) + voxels;
                dist = norm( (a-p) - ((a-p)*n')*n);
                
                if(dist < best_dist)
                    best_dist = dist;
                end
                
                % If the distance is larger than some value, we are outside
                % the channel.
                if(best_dist > center*0.2)
                    A(i,j,k) = 1;
                else
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
                    x(count) = i-center;
                    y(count) = j-center;
                    z(count) = k-center;
                    count = count + 1;
                end
            end
        end
    end

    %plot3(x,y,z,'ko')
    plot(x,y,'ko')
    axis('equal')
    
    fid = fopen('/projects/master/code/worlds/box_with_rotation_m.bin', 'w');
    fwrite(fid, N, 'uint');
    fwrite(fid, A, 'unsigned char');
    fclose(fid);
end