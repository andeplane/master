function testPlane

l = [ 1 3 10];
l0 = [-1 -1 -1];
line = [l; l0]

point = [1,2,3];
normal = [1,1,2];

%# a plane is a*x+b*y+c*z+d=0
%# [a,b,c] is the normal. Thus, we have to calculate
%# d and we're set
d = -point*normal'; %'# dot product for less typing

%# create x,y
[xx,yy]=ndgrid(-10:10,-10:10);

%# calculate corresponding z
z = (-normal(1)*xx - normal(2)*yy - d)/normal(3);

%# plot the surface
figure
surf(xx,yy,z)
hold on
plot3(line(:,1),line(:,2),line(:,3));
d = (point - l0)*normal'/(l*normal');

intersec = d*l + l0;
plot3(intersec(:,1),intersec(:,2),intersec(:,3),'ro');

end