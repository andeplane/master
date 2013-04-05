function plotting(vel,temp,pressure)
close all;

if temp
    temperature = dlmread('temperature.dat');

    figure(1)
    subplot(3,1,1);
    plot(temperature(:,1),temperature(:,2));
    legend('Temperature');
    xlabel('Time')
    ylabel('Temperature')
end

if vel
    % Open file
    velocity = dlmread('velocity.dat');

    % Get array size
    vel_size = size(velocity);
    % The last measured velocity
    v_x_last = velocity(vel_size(1),:);
    % Average across the canal
    avg_v = sum(v_x_last)/length(v_x_last)
    % Normalize
    %v_x_last = v_x_last/avg_v;
    
    subplot(3,1,2);
    plot(v_x_last);
    legend('Velocity profile');
    xlabel('x')
    ylabel('v')
end

if pressure
    pressure = dlmread('pressure.dat');

    subplot(3,1,3);
    x = 0:2;
    y = 0:1;
    imagesc(x,y,pressure);
    colorbar;

    % caxis([0 max(pressure)])

    legend('Pressure field');
    xlabel('x')
    ylabel('y')

    figure()

    world_bmp = imread('world_cool.bmp');
    world_bmp = sum(255-world_bmp,3); %Invert colors

    world_bmp(174,158,1) = 77;
    %world_bmp(174,158,2) = 0;
    %world_bmp(174,158,3) = 0;

    v_x = dlmread('vel_x.dat');
    v_y = dlmread('vel_y.dat');


    img_size = size(world_bmp)
    vel_size = size(v_x);

    img_size_x = img_size(2);
    img_size_y = img_size(1);

    vel_size_x = vel_size(2);
    vel_size_y = vel_size(1);

    x = linspace(1,img_size_x,vel_size_x);
    y = linspace(img_size_y,1,vel_size_y);
    mesh = meshgrid(x,y);

    imshow(world_bmp);
    hold on;
    quiver(x,y,v_x,-v_y,2);

    axis('equal');
end    
end