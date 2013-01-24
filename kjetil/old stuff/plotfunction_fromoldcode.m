%% Create movie:
close all
clear all
tic
foldername = 'Linear_channel_750_incl/';%Include /
%load ([foldername 'solution.mat']);
load([foldername 'output_iteration_number']);
tmax = output_iteration_number;
scrsz = get(0,'ScreenSize');
%h=figure('Position',[1 scrsz(4)/1.8 scrsz(3)/1.2 scrsz(4)/1.8])
%set(gcf,'renderer','zbuffer')
stepsize = 350;

for tind= stepsize:stepsize:tmax
    load([foldername 'data/' num2str(tind) '.mat']);
    display(['Step ' num2str(tind) ' of ' num2str(tmax)]);
    %figure(1)
    %clf
%     subplot(211)
%     gx = reshape(gcoord(1,Elem2Node(1:3,:)),3,nel);
%     gy = reshape(gcoord(2,Elem2Node(1:3,:)),3,nel);
%     hold off
%     patch(gx,gy,p)
% %     axis equal
% %     axis([-1 12 -1 1])
%     shading  interp
%     set(gcf,'renderer','zbuffer')
%     colorbar
%     colormap jet(500)
%     %set(gca,'CLim',[0 1500])
%     %title(['t = ' num2str(t(1)) ])
%     xlabel('x')
%     ylabel('y')
%     hold on
%     for i = 1:size(inclusions,2)
%         r = inclusions(3,i);
%         xc = inclusions(1,i);
%         yc = inclusions(2,i);
%         fnplt(rsmak('circle',r,[x y]),1,'w')
%         plot([-r*cos(inclusions(7,i)) r*cos(inclusions(7,i))]+x, [-r*sin(inclusions(7,i)) r*sin(inclusions(7,i))]+y,'w','linewidth',1);
%         
%         
% %             rotangle = inclusions(7,i);
% %             R = [cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
% %             phi = linspace(0,2*pi,50);
% %             x = r*cos(phi);
% %             y = r*sin(phi);
% %             rotpos = R*([x;y]);
% %             plot(rotpos(1,:)+xc,rotpos(2,:)+yc,'k','linewidth',3)
%         
%         
%     end
    
    
    
%     inclusions = inclusions(:,inclusions(1,:)>1.5);
%      inclusions = inclusions(:,inclusions(1,:)<5.5);
%     %subplot(212)
%     trisurf(Elem2Node(1:3,:)',gcoord(1,1:maxNode),gcoord(2,1:maxNode), u(1:ndim:ndim*maxNode))
%     axis equal
%     axis([2 5 -.35 .35])
% %     axis([9 10 -.25 .25])
%     grid off
%     colorbar
%     set(gca,'CLim',[0 6])
%     set(gcf,'visible','off');
%     colormap jet(500)
%     shading interp;
%     set(gcf,'renderer', 'zbuffer')
%     view(2)
%     xlabel('x')
%     ylabel('y')
%     hold on
%         for i = 1:size(inclusions,2)
%             r = inclusions(3,i);
%             x = inclusions(1,i);
%             y = inclusions(2,i);
%             theta = linspace(0,2*pi,20);
%             plot3(r*cos(theta)+x,r*sin(theta)+y,20*ones(length(theta)),'-k','linewidth',1);
%             %fnplt(rsmak('circle',r,[x y 1]),5,'w')
%             plot3([-r*cos(inclusions(7,i)) r*cos(inclusions(7,i))]+x, [-r*sin(inclusions(7,i)) r*sin(inclusions(7,i))]+y,[20 20],'k','linewidth',1);
%         
%         
%         
% 
% %             r = inclusions(3,i);
% %             xc = inclusions(1,i);
% %             yc = inclusions(2,i);
% %             theta = linspace(0,2*pi,20);
% %             rotangle = inclusions(7,i);
% %             R = [cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
% %             phi = linspace(0,2*pi,50);
% %             x = r*cos(phi);
% %             y = r*sin(phi)/2;
% %             rotpos = R*([x;y]);
% %             plot3(rotpos(1,:)+xc,rotpos(2,:)+yc,10*ones(size(x)),'k','linewidth',1)
% %         
%         
%         
%         
%         end
%         
        
        
        
        
        
        
        
    
    %sdrawnow
    %SaveMovie(:,tind)=getframe(h);
    %set(h, 'PaperPosition', [0 0 16 3.5]);
    %set(gcf, 'Papersize',[6 4]);
    %print(gcf,[foldername 'zoomed_figures/' num2str(round(tind/stepsize))],'-dpng','-r300')
    
    
    
    
%     
    nip = 6;
    [ipuv, ipw] = ip_triangle(nip);
    [N_arr, dNdu_arr] = shape_func_triangles(ipuv(:,1),ipuv(:,2),'tri7');
    [Pi_arr, dPidu_arr] = shape_func_triangles(ipuv(:,1),ipuv(:,2),'tri3');
    pressure = zeros(1,nel);
    AreaInc = zeros(1,nel);
    pressurecoord = zeros(2,nel);
    for ip = 1:nip
        dNdu = squeeze(dNdu_arr(ip,:,:));
        dPidu = squeeze(dPidu_arr(ip,:,:));
        Pi = Pi_arr(ip,:);
        N = N_arr(ip,:);
        for iel =1:nel
            index=Elem2Node(:,iel);
            xyel = gcoord(:,index);
            J = xyel*dNdu;
            detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
            AreaInc(iel) = AreaInc(iel) + ipw(ip)*detJ;
            
            index=Elem2Node(1:3,iel);
            xyel = gcoord(:,index);
            J = xyel*dPidu;
            detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
            pressure(iel) = pressure(iel) + ipw(ip)*detJ*Pi*p(:,iel);
        end
    end
    pressureInt(round(tind/stepsize)) = sum(pressure(gx(1,:)<2 & gx(1,:)>1.95))/sum(AreaInc( gx(1,:)<2 & gx(1,:)>1.95)) - sum(pressure(gx(1,:)>11))/sum(AreaInc(gx(1,:)>11))
    timeVal(round(tind/stepsize)) = t(1);
    
    
end



%save([foldername 'pressureInt.mat'],'pressureInt','timeVal');
toc
movie2avi(SaveMovie, [foldername 'movie.avi'])







%% Recalculate:


close all
clear all
plotvalues = {'.-b','.-r','.-m','.-c','.-k','.-g'}
plotname = {};
j=0;

dy = .1
edges = -1:dy:1;
histogramSumx = zeros(1,length(edges)-1);
histogramSumy = zeros(1,length(edges)-1);
startind = 1;
stopind = 781;

lowerx = -2;
upperx = 20;

for tind = startind:100:stopind
j=j+1;
foldername = 'Bagnolds_test5/';%Include /
load([foldername 'data/' num2str(tind) '.mat']);
load([foldername 'parameters.mat']);
saveflag = 0;


xpos = zeros(1,3*size(inclusions,2));

xpos(1:3:end) = inclusions(1,:);
xpos(2:3:end) = inclusions(2,:);
xpos(3:3:end) = inclusions(7,:);

%particleVel = EvaluateFunction(tind,xpos,inclusions,Mu,L_channel,h_channel,sineampl,channeltype,maxArea,np,reorder,method,gravity,saveflag,foldername);

x = inclusions(1,:);
y = inclusions(2,:);
vx = velx(2:end)./AreaInc(2:end);
vy = vely(2:end)./AreaInc(2:end);
vrot = curl./(2*AreaInc(1:end));

y = y(x>lowerx & x < upperx);
vx = vx(x>lowerx & x < upperx);
vy = vy(x>lowerx & x < upperx);
vrot = vrot(x>lowerx & x < upperx);
x = x(x>lowerx & x < upperx);


histogramy = zeros(1,length(edges)-1);
histogramx = zeros(1,length(edges)-1);

test = 0;
for i = 2:length(edges)
    test = test + length(vy (y<=edges(i) & y>edges(i-1) ));
    histogramy(i-1) = mean(vy (y<=edges(i) & y>edges(i-1) ));
    histogramx(i-1) = mean(vx (y<=edges(i) & y>edges(i-1) ));
end
figure(1)
hold on
plot(edges(2:end)-dy/2,histogramy,'Color',[0.9 0.9 0.9])
% plot(y,vy,'.m')
plotname{j} = ['t = ' num2str(t)];


figure(2)
hold on
plot(edges(2:end)-dy/2,histogramx,'Color',[0.9 0.9 0.9])


% figure(3)
% hold on
% plot(t,mean(vy),'.r')


% figure(4)
% % hold on
% colorvalue = jet(500);
% maxval = 1e-3;
% hold off
% for i = 1:length(x)
% plot(x(i),y(i),'.','Color',[colorvalue(round(vy(i)/maxval*250+250),:)],'MarkerSize',5)
% hold on
% end


drawnow
histogramSumy = histogramSumy+histogramy;
histogramSumx = histogramSumx+histogramx;
end

% legend(plotname);

% %% Plot
% figure(1)
% hold on
% edges = -.3:.025:.3;
% histogram = histc(y(x>6 & x < 11),edges);
% plot(edges,histogram/length(y),'r');


figure(1)
hold on
plot(edges(2:end)-dy/2,histogramSumy/(stopind-startind+1),'.-b','LineWidth',2)


figure(2)
hold on
plot(edges(2:end)-dy/2,histogramSumx/(stopind-startind+1),'.-b','LineWidth',2)

figure(4)
colorbar
set(gca,'clim',[-maxval,maxval])



%% Movement in y-direction
clear all
% close all

startind = 1;
stopind = 580;


j=1;
dx = 4;

for lowerx = 0:0
% lowerx = 0;
upperx = lowerx+dx;
% upperx=100;

foldername = 'Bagnold_test6/';%Include /
load([foldername 'data/' num2str(startind) '.mat']);
load([foldername 'parameters.mat']);

x = inclusions(1,:);
y = inclusions(2,:);
% vx = velx(2:end)./AreaInc(2:end);
% vy = vely(2:end)./AreaInc(2:end);


% vrot = curl;
plotindex = (x>lowerx & x < upperx);
y = y(plotindex);
% vx = vx(plotindex);
% vy = vy(plotindex);
% vrot = vrot(plotindex);
x = x(plotindex);


for tind = startind:stopind-startind-1:stopind
load([foldername 'data/' num2str(tind) '.mat']);
load([foldername 'parameters.mat']);

x2 = inclusions(1,:);
y2 = inclusions(2,:);


% vx2 = velx(2:end)./AreaInc(2:end);
% vy2 = vely(2:end)./AreaInc(2:end);
% vrot2 = curl./(2*AreaInc(1:end));
y2 = y2(plotindex);
% vx2 = vx2(plotindex);
% vy2 = vy2(plotindex);
% vrot2 = vrot2(plotindex);
x2 = x2(plotindex);


size(inclusions);


figure(1)
% hold on
plot(y,y2-y,'.')


% 
% figure(2)
% plot(x,y,'.')
% hold on
% plot(x2,y2,'.r')
% hold on

% plot([x;x2],[y;y2],'--k')



end

allcolors = lines(50);
figure(2)
hold on
plot(y,y2-y,'.','color','r','markersize',15)
j=j+1;

end


figure
plot(x,y,'.')
hold on
plot(x2,y2,'.r')
hold on

% plot([x;x2],[y;y2],'--k')



%% Analytical expression:
clear all
load flow_profile
vx = vx*0+2/3;
clear i
a = -1;

    
dy = y(2)-y(1);

b=0;



for j = 1:length(vx)-2
    

b = b   + (vx(j)*dy + (vx(j+1)-vx(j))*dy/2)


var1 = (((27*a^3 - 81*a - 81*b)^2 - 2916)^(1/2) + 27*a^3 - 81*a - 81*b)^(1/3);
var2 = 3*2^(1/3);


dx = [(var1/var2 + var2/var1 - a);...
    (-(1-i*sqrt(3))*var1/(2*var2) - 3*(1+i*sqrt(3))/(2^(2/3)*var1)-a); ...
    (-(1+i*sqrt(3))*var1/(2*var2) - 3*(1-i*sqrt(3))/(2^(2/3)*var1)-a)]


% var1 = (((-27*a^3 + 81*a + 81*b)^2 - 2916)^(1/2)-27*a^3 + 81*a + 81*b)^(1/3);
% var2 = 3*2^(1/3);
% dx = [(-var1/var2 - var2/var1 - a);...
%     ((1-i*sqrt(3))*var1/(2*var2) + 3*(1+i*sqrt(3))/(2^(2/3)*var1)-a); ...
%     ((1+i*sqrt(3))*var1/(2*var2) + 3*(1-i*sqrt(3))/(2^(2/3)*var1)-a)];





% displ(j) = real(dx(real(dx)<=2 & real(dx)>=0 & imag(dx)<1e-3))

displ(j) = real(dx(3));

displ0(j) = dy*j;


end

hold on
plot(y(2:end-1),(displ-displ0),'.-k')









%% Volume fraction in channel:


stepsize = 3;
tmax = 500;
h_channel = .5;
ampl = .1;
volumeFrac = zeros(1,round(tmax/stepsize)-1);
foldername = 'Linear_channel_750_incl/';%Include /
for tind= stepsize:stepsize:tmax
    load([foldername 'data/' num2str(tind) '.mat']);
    inclusions = inclusions(:,inclusions(1,:)>2 & inclusions(1,:)<11);
    if size(inclusions,2)>0
        volumeFrac(round(tind/stepsize)) = sum(pi*inclusions(3,:).^2);
    end
    size(inclusions,2)
end


% x_rough = linspace(0,9,length(rough_surface));
% load surf_rough_top;
% total_rough = rough_surface*ampl;
% load surf_rough_bottom;
% total_rough = total_rough - rough_surface*ampl;
% dx = x_rough(2)-x_rough(1);
% totalvolume = sum(total_rough.*dx) + h_channel*9


totalvolume = h_channel*9;
volumeFrac = volumeFrac/totalvolume;
save([foldername 'volumeFrac.mat'],'volumeFrac');
figure
plot(volumeFrac)


%% Particle curl:


clear all
close all
clc
channelstart=2;
channelend=10;
tmax = 380;
curl_tot = zeros(1,tmax);
curl_max = curl_tot;
timeVal = curl_tot;
edges = -.2:.005:.2;
histogram_curl = zeros(size(edges));
foldername = 'Linear_channel_100_incl_2/';%Include /
for tind= 100:tmax
    display(num2str(tind))
    load([foldername 'data/' num2str(tind) '.mat']);
    %curl=abs(curl(2:end))./2;
    curl = curl./2;
    curl = curl(inclusions(1,:)>channelstart & inclusions(1,:)<channelend);
    if length(curl)>0
        curl_max(tind)=max(curl);
        curl_tot(tind)=mean(curl);
    end
    timeVal(tind)=t(1);
    histogram_curl = histogram_curl + histc(curl,edges);
end


figure(1)
subplot(211)
plot(timeVal,curl_tot)

subplot(212)
plot(timeVal,curl_max)


figure(2)
plot(edges,histogram_curl/sum(histogram_curl))






figure(3)
subplot(3,1,1:2)
hold on
plot(edges,histogram_curl/sum(histogram_curl),'b')
axis([-.1 .1 0 .3])
xlabel('$\omega$')
ylabel('count / total count')
legend('Circle','Ellipse')

subplot(313)
hold on
plot(timeVal,curl_max,'b')
xlabel('t')
ylabel('$max(|\omega|)$')
axis([0 max(t) 0 .15])



histogram_curl = zeros(size(edges));
curl_tot = zeros(1,tmax);
foldername = 'Linear_channel_ellipse_100_incl/';%Include /
for tind= 100:tmax
    display(num2str(tind))
    load([foldername 'data/' num2str(tind) '.mat']);
    %curl=abs(curl(2:end))./2;
    curl = curl./2;
    curl = curl(inclusions(1,:)>channelstart & inclusions(1,:)<channelend);
    if length(curl)>0
        curl_max(tind)=max(curl);
        curl_tot(tind)=mean(curl);
    end
    timeVal(tind)=t(1);
    histogram_curl = histogram_curl + histc(curl,edges);
end


figure(1)
subplot(211)
hold on
plot(timeVal,curl_tot,'r')
xlabel('t')
ylabel('$<|\omega|>$')
axis([0 max(t) 0 .04])

subplot(212)
hold on
plot(timeVal,curl_max,'r')
xlabel('t')
ylabel('$max(|\omega|)$')
axis([0 max(t) 0 .15])


figure(2)
hold on
plot(edges,histogram_curl/sum(histogram_curl),'r')
axis([-.1 .1 0 .3])
xlabel('$\omega$')
ylabel('count / total count')
legend('Circle','Ellipse')




figure(3)
subplot(3,1,1:2)
hold on
plot(edges,histogram_curl/sum(histogram_curl),'r')
axis([-.1 .1 0 .3])
xlabel('$\omega$')
ylabel('count / total count')
legend('Circle','Ellipse')

subplot(313)
hold on
plot(timeVal,curl_max,'r')
xlabel('t')
ylabel('$max(|\omega|)$')
axis([0 max(t) 0 .15])

%%

foldername = 'Linear_channel_ellipse_100_incl/';%Include /
load([foldername 'data/' num2str(200) '.mat']);
totalrot = inclusions(7,:);
load([foldername 'data/' num2str(tmax) '.mat']);
totalrot = totalrot - inclusions(7,:);

%% Value in set of points:
% clear all; close all; clc;

intN = 1;


% 
% Npoints = 1000;
% points_all = [rand(Npoints,1)*10+5 rand(Npoints,1)*2-1]

 [X,Y]=meshgrid(linspace(-1,0,10),linspace(-1,1,100));
 

points_all = [X(:) Y(:)];
Npoints = length(X(:));

vx_local = zeros(1,Npoints);
vy_local = zeros(1,Npoints);

foldername = 'Bagnolds_test_ordered3/';%Include /
load([foldername 'data/' num2str(intN) '.mat']);
load([foldername 'parameters.mat']);
ux = u(1:2:end);
uy = u(2:2:end);

for i = 1:size(points_all,1)
    
display([num2str(i) '/' num2str(Npoints)])
points = points_all(i,:);

%Find element
TRIindex = Elem2Node(1:3,:); TRIindex=unique(TRIindex(:));
t = tsearchn(gcoord(:,TRIindex)',Elem2Node(1:3,:)',points);
ux_local = ux(Elem2Node(:,t));
uy_local = uy(Elem2Node(:,t));

%Transfer to local coordinates:
xyel = gcoord(:,Elem2Node(:,t));
points_local = [points(:,1)' - xyel(1,1); ...
                points(:,2)' - xyel(2,1)];
xyel = [xyel(1,:) - xyel(1,1); ...
       xyel(2,:) - xyel(2,1)];
A = [xyel(1,2) xyel(1,3); ...
     xyel(2,2) xyel(2,3)];
uvel = inv(A)*xyel;
points_local = inv(A)*points_local;

%Calculate exact value from shape functions
[N, ~]    = shp_deriv_triangle(points_local', 7);
vx_local(i) = ux_local'*N{1};
vy_local(i) = uy_local'*N{1};

end

    

for j = 1:size(Y,1)
vx_sum(j) = sum(vx_local(points_all(:,2)==Y(j,1)))/size(Y,2)
vy_sum(j) = sum(vy_local(points_all(:,2)==Y(j,1)))/size(Y,2)
end

figure(1)
% plot(points_all(:,2),vx_local,'.')
hold on
plot(Y(:,1),vx_sum,'-r')
% plot(Y(:,1),vy_sum,'-r')

x = linspace(-1,1,1000);
hold on
plot(x,1-x.^2,'--k')

y = Y(:,1);
vx = vx_sum;
save('flow_profile.mat','vx','y');

sum(vx)*(y(2)-y(1))


