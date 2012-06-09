function plotEnergy(t,E)
    figure(3)
    plot(t*10^9,E)
    xlabel('t [ns]');
    ylabel('E [J]');
    title('Energy');
end