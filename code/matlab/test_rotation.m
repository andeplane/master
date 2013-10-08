% Creates a box with a slice
function test_rotation(theta, voxels, i, j)
    figure;
    
    create_box_with_rotation(theta, voxels, 'ko');
    %hold on
    %create_box_with_rotation(0, voxels, 'r.');
    center = round(voxels / 2);
    sprintf('Center at : %f', center)
    sprintf('dy limit: %f', center*0.2)
    
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    dx = i - center;
    dy = j - center;
    v = [dx dy]'
    %plot(v(1), v(2),'ro');
    v = R*v
    %plot(v(1), v(2),'bo');
    if v(2) < -center
        disp('Too negative, increasing...')
        v(2) = v(2) + center;
    end

    if v(2) > center
        disp('Too positive, decreasing...')
        v(2) = v(2) - center;
    end
    v
    if(abs(v(2)) > center*0.2)
        disp('Voxels is filled')
    else
        disp('Voxels is not filled')
    end
    
    