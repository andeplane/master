function animate(numParticles,p,r,R,L,i,tau,fig1,circleX,circleY,winsize,a)
    % Clear figure and plot this timestep
    particleIndices = p(1:numParticles);
    figure(1)
    subplot(2,1,1)
    hold off
    plot(r(particleIndices,1),r(particleIndices,2),'d');
    hold on;
    plot([0,L],[-R,-R],'g')
    plot([0,L],[R,R],'g')
    % plot(circleX,circleY,'g');
    xlabel('z [m]');
    ylabel('y [m]');
    axis('equal')

    %axis([-1.1*R 1.1*R -1.1*R 1.1*R]);
    axis([0 L -1.1*R 1.1*R]);

    title(sprintf('t = %f ns',10^9*i*tau));
    
    subplot(2,1,2)
    
    hold off
    plot(r(particleIndices,3),r(particleIndices,2),'d');
    hold on;
    
    plot(circleX,circleY,'g');
    xlabel('z [m]');
    ylabel('y [m]');
    axis('equal')

    % axis([-1.1*R 1.1*R -1.1*R 1.1*R]);
    
    title(sprintf('t = %f ns',10^9*i*tau));    
end