function createMovie()

%% Create movie:
close all
tic
foldername = 'Bagnold_test6/';%Include /
%load ([foldername 'solution.mat']);
load([foldername 'output_iteration_number']);
tmax = output_iteration_number-1;
scrsz = get(0,'ScreenSize');
h=figure('Position',[1 scrsz(4)/2.2 scrsz(3)/1 scrsz(4)/2.2])

stepsize = 1;
plotinc = true;
plotfield = true
savemovie = false;
printpng = false;

myVideo = VideoWriter([foldername 'myfile.avi']);
% myVideo.FrameRate = 30;  % Default 30
% myVideo.Quality = 50;    % Default 75
open(myVideo)
j=1;

for tind= 1:15:tmax+1%stepsize:stepsize:tmax
        display(['Step ' num2str(tind) ' of ' num2str(tmax)]);

    if plotfield
    load([foldername 'data/' num2str(tind) '.mat']);
    gx = reshape(gcoord(1,Elem2Node(1:3,:)),3,nel);
    gy = reshape(gcoord(2,Elem2Node(1:3,:)),3,nel);
    figure(1)
    hold off
    trisurf(Elem2Node(1:3,:)',gcoord(1,1:maxNode),gcoord(2,1:maxNode), u(1:ndim:ndim*maxNode))
    axis equal
    grid off
    shading interp;
    set(gcf,'renderer', 'zbuffer')
    view(2)
%     set(gca,'clim',[0 .15])
%     set(gca,'clim',[0 .125])
% set(gca,'clim',[0 1.1])
    colorbar
    hold on
    end
    
    if plotinc
    for i = 1:size(inclusions,2)
        r = inclusions(3,i);
        x = inclusions(1,i);
        y = inclusions(2,i);
        theta = linspace(0,2*pi,20);
        plot3(r*cos(theta)+x,r*sin(theta)+y,10*ones(length(theta)),'-k','linewidth',.8);
        %fnplt(rsmak('circle',r,[x y 1]),5,'w')
        plot3([-r*cos(inclusions(7,i)) r*cos(inclusions(7,i))]+x, [-r*sin(inclusions(7,i)) r*sin(inclusions(7,i))]+y,[10 10],'k','linewidth',.8);
    end
%     axis([2 20 -1.2 1.2])
%     axis off

    end
    drawnow
    
%     
%     figure(2)
%     patch(gx,gy,p)
%     axis equal
%     shading interp
%     set(gcf,'renderer','zbuffer')
% %     set(gca,'clim',[0 1e9])
% %     axis([3 10 -.5 .5])
%     colorbar
    

    drawnow

    if printpng
    set(h, 'PaperPosition', [1 scrsz(4)/4 scrsz(3)/1 scrsz(4)/4]./100);
    print(h,[foldername 'figures/' num2str(round(tind/stepsize))],'-dpng','-zbuffer','-r400')
    end
    
    
    
    if savemovie
    img = imread([foldername 'figures/' num2str(round(tind/stepsize)) '.png']);
    writeVideo(myVideo,img)
    end
    

end


close(myVideo)


% movie2avi(mov, 'myPeaks.avi', 'compression', 'None');

toc