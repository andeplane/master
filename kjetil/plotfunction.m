
clc
clear all;
load mesh.mat;
close all
fclose all;

for i = 0:1:10000
   
[~,~,~,u] = read_vtk(sprintf('output/v%06d.vtu',i));
[~,~,~,p] = read_vtk(['output/p' num2str(i,'%06d') '.vtu']);
vx = u(1:3:end);
vy = u(2:3:end);
vz = u(3:3:end);

    %Plot~
    figure(1)
    hold off
    % plot3(gcoord(1,:),gcoord(2,:),vy,'.')
    trisurf(Elem2Node',gcoord(1,:),gcoord(2,:),p)
    axis equal
    grid off
    shading interp;
    set(gcf,'renderer', 'zbuffer')
    view(2)
    %set(gca,'clim',[0 0.15])
    colorbar
    hold on
    quiver3(gcoord(1,:),gcoord(2,:),gcoord(1,:)*0+10, vx, vy, vy*0,'k');
    hold off
    drawnow
end