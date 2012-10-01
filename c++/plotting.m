function plotting
clear all;
close all;

%energy = dlmread('energy.dat');
temperature = dlmread('temperature.dat');

%t = energy(:,1);
%Ek = energy(:,2);
%Ep = energy(:,3);
%E = energy(:,4);

figure(1)
%subplot(2,1,1);
%plot(t,Ek,'r');
%hold on
%plot(t,Ep,'g');
%plot(t,E,'b');
%legend('Kinetic energy','Potential energy','Total energy');
%xlabel('Time')
%ylabel('Energy')

plot(temperature(:,1),temperature(:,2));
legend('Temperature');
xlabel('Time')
ylabel('Temperature')
hold on
x = zeros(length(temperature(:,1)))+3;
plot(x,'r');
    
end