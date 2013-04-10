function create_cylinder(fraction)

voxels = 100;
N = [voxels voxels voxels];
A = zeros(N);
center = 50;

cylinder_radius = center*fraction;

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

x = zeros(1000,1);
y = zeros(1000,1);
z = zeros(1000,1);
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
fid = fopen('/projects/master/code/worlds/cylinder_m.bin', 'w');
fwrite(fid, N, 'unsigned char');
fwrite(fid, A, 'unsigned char');
fclose(fid);