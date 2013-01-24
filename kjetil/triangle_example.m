clear all; close all;

coord = [-1 -1 1  1;...
        -1  1 1 -1];




segm = [2 4;...
    3 1;...
    2 3]; % connect nodes
segm = [segm [1 3; 2 4; 0 0]];

hole = [];
regions = [-1+1e-3;...
    -1+1e-3;...
    1];

%Circle:
n = 8;
phi = linspace(0,2*pi,n+1);
phi(end) = [];
xc = 0;
yc = 0;
r = 0.2;
x = r*cos(phi)+xc;
y = r*sin(phi)+yc;
segmxy = [size(segm,2)+1:size(segm,2)+n; ...
    size(segm,2)+2:size(segm,2)+n+1; ...
    5+zeros(1,n)];
coord   = [coord [x;y]];
segm    = [segm segmxy];
segm(2,end) = size(segm,2)+1-n;
hole = [xc+r-1e-1; yc; 10];





figure
plot(coord(1,:),coord(2,:),'.')


[gcoord,Elem2Node,bc_node,phase_id] = triangle_mesh_generation(coord,segm,hole,regions,6);



%% Periodic boundaries:



xr = gcoord(1,bc_node(1,bc_node(2,:)==1));
yr = gcoord(2,bc_node(1,bc_node(2,:)==1));
yr = yr(xr>0);
xr = xr(xr>0);
yl = yr;
xl = xr-2;
[yl I] = sort(yl);
xl = xl(I);
[yr I] = sort(yr,'descend');
xr = xr(I);
% Create new mesh:
coord = [-1 xl -1 1 xr 1;...
        -1 yl 1 1 yr -1];

n = length(coord);
segm = [1:n; ...
    2:n 1; ...
    zeros(1,length(xl)+1) 2 zeros(1,length(xl)+1) 3];
figure
plot(coord(1,:),coord(2,:),'.-')


hole = [];
regions = [1-1e-3;1-1e-3;1];




% %Circle:
% n = 8;
% phi = linspace(0,2*pi,n+1);
% phi(end) = [];
% xc = 0;
% yc = 0;
% r = 0.2;
% x = r*cos(phi)+xc;
% y = r*sin(phi)+yc;
% segmxy = [size(segm,2)+1:size(segm,2)+n; ...
%     size(segm,2)+2:size(segm,2)+n+1; ...
%     5+zeros(1,n)];
% coord   = [coord [x;y]];
% segm    = [segm segmxy];
% segm(2,end) = size(segm,2)+1-n;
% hole = [xc+r-1e-1; yc; 10];


[gcoord,Elem2Node,bc_node,phase_id] = triangle_mesh_generation(coord,segm,hole,regions,6);







%% Plot
figure;
clf
plot(gcoord(1,:),gcoord(2,:),'.b');
hold on;
nel=size(Elem2Node,2);
gx = reshape(gcoord(1,Elem2Node(1:3,:)),3, nel);
gy = reshape(gcoord(2,Elem2Node(1:3,:)),3, nel);

patch(gx(:,phase_id==1),gy(:,phase_id==1),'r')
patch(gx(:,phase_id==2),gy(:,phase_id==2),'y')
axis equal;



% 
% 
% bc_node(2,1)=2;
% bc_node(2,2)=2;
% bc_node(2,3)=2;
% bc_node(2,4)=2;

figure
plotsegm=3;
plot(gcoord(1,bc_node(1,bc_node(2,:)==plotsegm)),gcoord(2,bc_node(1,bc_node(2,:)==plotsegm)),'.r')
axis equal
save('Meshdata', 'gcoord', 'Elem2Node' ,'bc_node','phase_id')

