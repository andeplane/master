function collide_with_voxel(r, v, Lx, Ly, Lz, i, j, k, Nx, Ny, Nz)
    format long;
    voxel_size = [Lx / Nx, Ly / Ny, Lz / Nz];
    normals = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
    
    p1 = [i*voxel_size(1), j*voxel_size(2), k*voxel_size(3)];
    p2 = [(i+1)*voxel_size(1), j*voxel_size(2), k*voxel_size(3)];
    p3 = [i*voxel_size(1), j*voxel_size(2), (k+1)*voxel_size(3)];
    p4 = [(i+1)*voxel_size(1), j*voxel_size(2), (k+1)*voxel_size(3)];
    p5 = [i*voxel_size(1), (j+1)*voxel_size(2), k*voxel_size(3)];
    p6 = [(i+1)*voxel_size(1), (j+1)*voxel_size(2), k*voxel_size(3)];
    p7 = [i*voxel_size(1), (j+1)*voxel_size(2), (k+1)*voxel_size(3)];
    p8 = [(i+1)*voxel_size(1), (j+1)*voxel_size(2), (k+1)*voxel_size(3)];
    
    plot3(p1(1), p1(2), p1(3),'ko');
    hold on
    plot3(p2(1), p2(2), p2(3),'ko');
    plot3(p3(1), p3(2), p3(3),'ko');
    plot3(p4(1), p4(2), p4(3),'ko');
    plot3(p5(1), p5(2), p5(3),'ko');
    plot3(p6(1), p6(2), p6(3),'ko');
    plot3(p7(1), p7(2), p7(3),'ko');
    plot3(p8(1), p8(2), p8(3),'ko');
    
    time_1 = time_until_collision(r, v, p1, normals(1,:) )
    time_2 = time_until_collision(r, v, p2, normals(2,:) )
    time_3 = time_until_collision(r, v, p1, normals(3,:) )
    time_4 = time_until_collision(r, v, p5, normals(4,:) )
    time_5 = time_until_collision(r, v, p1, normals(5,:) )
    time_6 = time_until_collision(r, v, p3, normals(6,:) )
    
    isp1 = r + v*time_1;
    isp2 = r + v*time_2;
    isp3 = r + v*time_3;
    isp4 = r + v*time_4;
    isp5 = r + v*time_5;
    isp6 = r + v*time_6;
    
    plot3([r(1) isp1(1)], [r(2) isp1(2)], [r(3) isp1(3)], 'ro-')
    plot3([r(1) isp2(1)], [r(2) isp2(2)], [r(3) isp2(3)], 'ro-')
    plot3([r(1) isp3(1)], [r(2) isp3(2)], [r(3) isp3(3)], 'ro-')
    plot3([r(1) isp4(1)], [r(2) isp4(2)], [r(3) isp4(3)], 'ro-')
    plot3([r(1) isp5(1)], [r(2) isp5(2)], [r(3) isp5(3)], 'ro-')
    plot3([r(1) isp6(1)], [r(2) isp6(2)], [r(3) isp6(3)], 'ro-')
    
    plot3(r(1), r(2), r(3), 'go');
    xlim([p1(1) p2(1)]);
    ylim([p1(1) p5(1)]);
    zlim([p1(1) p3(1)]);

end

function dot_prod = time_until_collision(r, v, point_in_plane, normal_vector)
    dot_prod = dot( (point_in_plane - r) , normal_vector) / dot(v,normal_vector);
end
