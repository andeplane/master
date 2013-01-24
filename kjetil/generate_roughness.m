


ampl = .1;
numb = 9;
n_peak = 5;
H = 0.9;

%meanvalue = 0;
meanmax = 5;
meanall = zeros(1,meanmax);
for meancounter = 1:meanmax
N = 2^numb+1;
a_0 = 1/N;
dy = 1;
mid = N;
surf = zeros(N,1);
closed = zeros(N,1);
dx = N-1;
dx = dx/2;
for i = 2:numb+1
    for j = dx+1:dx:N-1
        if(closed(j) ~= 1)
            if i <= n_peak
                losed(j) = 1;
            end
            if i > n_peak
            closed(j) = 1;
            surf(j) = (surf(j-dx) + surf(j+dx))/(2);
            surf(j) = (surf(j) + (rand-0.5)*(a_0*dx)^(H));
            end
        end
    end
    dx = dx/2;
end
end

rough_surface=surf/max(surf);
plot(linspace(0,1,N),surf)
title(['$H = $' num2str(H)])



save('surf_rough_bottom.mat','rough_surface')


writeMeshToFile
