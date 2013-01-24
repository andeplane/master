%function [gcoord, Elem2Node ,bc_node,phase_id,inclusions] = createmesh(inclusions,maxArea,L,h,type,sineampl,foldername,particletype)
%inclusions(N_incl,4): xc, yc, r, n, phase_id
clc
clear all
close all

maxArea = 0.1;
L = 10;
h = 1;
foldername = '';
type = 'rough';
sineampl = .25;

hole = [];
regions = [0;...
    -1+1e-5;...
    1];

%Channel properties:
left_box = [1 -2 -2 1;-1 -1 1 1];
switch type
    case 'linear'
        narrowchannel_top =[2  1+L ;
            h/2 h/2];
        narrowchannel_bottom =[1+L   2;
            -h/2 -h/2];

    case 'linear2'
        narrowchannel_top =[2;
            h/2];
        narrowchannel_bottom =[2;
            -h/2];
        
    case 'linear_decreasing'
        narrowchannel_top =[2  1+L ;
            h/2 h/2-sineampl/2];
        narrowchannel_bottom =[1+L   2;
            -h/2+sineampl/2 -h/2];
        
    case 'sine'
        narrowchannel_top =[linspace(2,1+L,100);
            sineampl*cos(linspace(0,L*2*pi,100))+h/2];
        narrowchannel_bottom =[linspace(1+L,2,100);
            sineampl*cos(linspace(L*2*pi,0,100))-h/2];
        
        
    case 'sine_opposite'
        narrowchannel_top =[linspace(2,1+L,100);
            sineampl*cos(linspace(0,L*2*pi,100))+h/2];
        narrowchannel_bottom =[linspace(1+L,2,100);
            -sineampl*cos(linspace(L*2*pi,0,100))-h/2];
        
    case 'rough'
        load surf_rough_top.mat
        narrowchannel_top =[linspace(2,1+L,length(rough_surface));
            sineampl*rough_surface'+h/2];
        load surf_rough_bottom.mat
        narrowchannel_bottom =[linspace(1+L,2,length(rough_surface));
            sineampl*rough_surface'-h/2];
end

right_box = [1+L 4+L 4+L 1+L; 1 1 -1 -1];

coord = [left_box narrowchannel_top right_box narrowchannel_bottom];

segm = [1:length(coord); ...
    2:length(coord) 1; ...
    1 2 3 4 4*ones(1,size(narrowchannel_top,2)) ...
    5 ...
    6 ...
    5 ...
    4*ones(1,size(narrowchannel_bottom,2)+1)];

% segm = [1:length(coord); ...
%     2:length(coord) 1; ...
%     1 2 3 4 4*ones(1,size(narrowchannel_top,2)) ...
%     5*ones(1,size(right_box,2)-1) ...
%     ...
%     ...
%     4*ones(1,size(narrowchannel_bottom,2)+1)];

%% Add inclusions:
if false
switch particletype
    
    case 'circle'
        for i = 1:length(inclusions(1,:))
            xc = inclusions(1,i);
            yc = inclusions(2,i);
            r = inclusions(3,i);
            n = inclusions(4,i);
            id = inclusions(5,i);
            phi = linspace(0,2*pi,n+1);
            phi(end) = [];
            x = r*cos(phi)+xc;
            y = r*sin(phi)+yc;
            segmxy = [size(segm,2)+1:size(segm,2)+n; ...
                size(segm,2)+2:size(segm,2)+n+1; ...
                10+zeros(1,n)];
            coord   = [coord [x;y]];
            segm    = [segm segmxy];
            segm(2,end) = size(segm,2)+1-n;
            regions = [regions [xc; yc; id]];
        end
        
    case 'ellipse'
        for i = 1:length(inclusions(1,:))
            
            rotangle = inclusions(7,i);
            R = [cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
            xc = inclusions(1,i);
            yc = inclusions(2,i);
            r = inclusions(3,i);
            n = inclusions(4,i);
            id = inclusions(5,i);
            phi = linspace(0,2*pi,n+1);
            phi(end) = [];
            x = r*cos(phi);
            y = r*sin(phi)/2;
            rotpos = R*([x;y]);
            x = rotpos(1,:)+xc;
            y = rotpos(2,:)+yc;
            segmxy = [size(segm,2)+1:size(segm,2)+n; ...
                size(segm,2)+2:size(segm,2)+n+1; ...
                10+zeros(1,n)];
            coord   = [coord [x;y]];
            segm    = [segm segmxy];
            segm(2,end) = size(segm,2)+1-n;
            regions = [regions [xc; yc; id]];
            
        end
        
        
        
    case 'car'
        for i = 1:length(inclusions(1,:))
            
            rotangle = inclusions(7,i);
            R = [cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
            xc = inclusions(1,i);
            yc = inclusions(2,i);
            r = inclusions(3,i);
            n = inclusions(4,i);
            id = inclusions(5,i);
            phi = linspace(0,2*pi,n+1);
            phi(end) = [];
            x = r*cos(phi);
            y = r*sin(phi)/2;
            rotpos = R*([x;y]);
            x = rotpos(1,:)+xc;
            y = rotpos(2,:)+yc;
            segmxy = [size(segm,2)+1:size(segm,2)+n; ...
                size(segm,2)+2:size(segm,2)+n+1; ...
                10+zeros(1,n)];
            coord   = [coord [x;y]];
            segm    = [segm segmxy];
            segm(2,end) = size(segm,2)+1-n;
            regions = [regions [xc; yc; id]];
            
            
            
        end
        
        
end
end

%% Generate mesh:
[gcoord,Elem2Node,bc_node,phase_id] = triangle_mesh_generation(coord,segm,hole,regions,3,'',maxArea,foldername);

%% Find boundary nodes:
%bc_b_lbox = bc_node(1,bc_node(2,:)==1);
%bc_l_lbox = bc_node(1,bc_node(2,:)==2);
%bc_t_lbox = bc_node(1,bc_node(2,:)==3);
%bc_channel = bc_node(1,bc_node(2,:)==4);
%bc_right_box = bc_node(1,bc_node(2,:)==5);

%% Plot
figure;
clf
nel=size(Elem2Node,2);
gx = reshape(gcoord(1,Elem2Node(1:3,:)),3, nel);
gy = reshape(gcoord(2,Elem2Node(1:3,:)),3, nel);
patch(gx(:,phase_id==1),gy(:,phase_id==1),'w')
patch(gx(:,phase_id>=2),gy(:,phase_id>=2),'b')
axis equal
set(gcf,'renderer','zbuffer')
grid off
drawnow

%%
writeMeshToFile
save('mesh.mat');

