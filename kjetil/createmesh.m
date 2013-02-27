%function [gcoord, Elem2Node ,bc_node,phase_id,inclusions] = createmesh(inclusions,maxArea,L,h,type,sineampl,foldername,particletype)
%inclusions(N_incl,4): xc, yc, r, n, phase_id
clc
clear all
close all

maxArea = 0.01;
L = 10;
h = 2;
foldername = '';
type = 'linear2';
sineampl = 0;
N = 500;

hole = [];
regions = [-2+1e-16;...
    -1+1e-16;...
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
        narrowchannel_top =[linspace(2,L+1,N);
            sineampl*cos(linspace(0,L*2*pi,N))+h/2];
        narrowchannel_bottom =[linspace(L+1,2,N);
            sineampl*cos(linspace(L*2*pi,0,N))-h/2];
        
        
    case 'sine_opposite'
        narrowchannel_top =[linspace(2,1+L,1000);
            sineampl*cos(linspace(0,L*2*pi,1000))+h/2];
        narrowchannel_bottom =[linspace(1+L,2,1000);
            -sineampl*cos(linspace(L*2*pi,0,1000))-h/2];
        
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
    3 1 3 3 3*ones(1,size(narrowchannel_top,2)) ...
    3 ...
    2 ...
    3 ...
    3*ones(1,size(narrowchannel_bottom,2)+1)];


%% Add inclusions:
if true
    inclusions = [0 0 .2 50 4]';
    particletype = 'circle';
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
                    id+zeros(1,n)];
                coord   = [coord [x;y]];
                segm    = [segm segmxy];
                segm(2,end) = size(segm,2)+1-n;
                hole = [hole [xc; yc;id]];
                %regions = [regions [xc; yc; id]];
            end
            
            
            
        case 'ellipse'
            display('fsa')
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

%% Find facets:
tic
[Elem2Node,I]=sort(Elem2Node,'ascend');
FacetID0 = [];
Elemind = 1:size(Elem2Node,2);
bc_test=[];
for ID = 1:max(bc_node(2,:))
bc_test = [bc_test bc_node(1,bc_node(2,:)==ID)];
end
Elemtest = zeros(size(Elem2Node));
for i = 1:length(bc_test)
    Elemtest = Elemtest + (Elem2Node==bc_test(i));
end
Cellindex = Elemind(sum(Elemtest,1)>1);
Facetindex = Elemtest(1,Cellindex)+2*Elemtest(2,Cellindex)+3*Elemtest(3,Cellindex)-2;

ind = 1:size(bc_node,2);
I2 = zeros(1,size(bc_node,2));
for i = 1:size(bc_node,2)
    a=Elem2Node==bc_node(1,i);
    b=Elemind(sum(a)>0);
    for j = 1:length(b)
       if sum(Cellindex==b(j))==1
        add(ind(Cellindex==b(j)))=bc_node(2,i);
       end
    end
end
%[a2,I]=sort(I2);
%add = bc_node(2,I);

FacetID0 = [Cellindex; Facetindex; add]';
FacetID0 = FacetID0(FacetID0(:,2)<4,:);
counter = repmat(Elemind,3,1);
counter2 = repmat(1:3,1,length(Elemind));
FacetID = [counter(:) counter2(:) zeros(size(Elem2Node,2)*3,1)];
FacetID(FacetID0(:,1)*3 - FacetID0(:,2)+1 ,3) = FacetID0(:,3);
FacetID(:,1:2)=FacetID(:,1:2)-1;
display(['Find facets: ' num2str(toc) 's'])




%% Plot
figure;
clf
nel=size(Elem2Node,2);
gx = reshape(gcoord(1,Elem2Node(1:3,:)),3, nel);
gy = reshape(gcoord(2,Elem2Node(1:3,:)),3, nel);
patch(gx(:,Cellindex),gy(:,Cellindex),'w')
%patch(gx(:,phase_id>=2),gy(:,phase_id>=2),'b')
axis equal
set(gcf,'renderer','zbuffer')
grid off

figure;
plotcolor = jet(max(bc_node(2,:)));
hold on
for i = 1:max(bc_node(2,:))
    plot(gcoord(1,bc_node(1,bc_node(2,:)==i)),gcoord(2,bc_node(1,bc_node(2,:)==i)),'.','color',plotcolor(i,:))
    patch(gx(:,FacetID0(:,1)),gy(:,FacetID0(:,1)),'w')
end
drawnow

%%
writeMeshToFile
save('mesh.mat');

