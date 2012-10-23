function plotting
clear all;
close all;

%energy = dlmread('energy.dat');
temperature = dlmread('temperature.dat');

velocity = dlmread('velocity.dat');

vel_size = size(velocity);
v_x_last = velocity(vel_size(1),:);
avg_v = sum(v_x_last)/length(v_x_last);
v_x_last = v_x_last/avg_v;

v_x_last = smooth(smooth(smooth(v_x_last)));


figure(1)
subplot(2,1,1);
plot(temperature(:,1),temperature(:,2));
legend('Temperature');
xlabel('Time')
ylabel('Temperature')

subplot(2,1,2);
plot(v_x_last);
legend('Velocity profile');
xlabel('x')
ylabel('v')
    
end