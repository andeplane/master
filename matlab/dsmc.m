function dsmc(timesteps)
    close all;
    format long;
    
    % Plots
    energyVsTime = true;
    
    particlesPerCell = 25;
    particles = 5000;
    
    % Calculate the required number of cells
    cells = ceil((particles/particlesPerCell)^(1/3));
    
    % Some constants
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    mass = 6.63e-26;       % Mass of argon atom (kg)
    diam = 3.66e-10;       % Effective diameter of argon atom (m)
    T = 273;               % Temperature (K)
    density = 1.78;        % Density of argon at STP (kg/m^3)
    L = 1e-5;              % Tube length
    R = L/10;              % Tube radius
    
    volume = pi*R*R*L;
    
    eff_num = density/mass*volume/particles;
    v0 = sqrt(3*boltz*T/mass); % Used to calculate initial vmax
    sigma = sqrt(boltz*T/mass); % Standard deviation
    
    % For small systems (testing) the number of cells^2 might be larger
    % than the number of particles
    dim = max(100*particles,cells^3) + 10;
    sortData = zeros(dim,3);
    
    %sortData explanation:
       %cell_n = sortData(:,1) Number of particles in a specific cell
       %index = sortData(:,2) Ordered set of particles
       %Xref = sortData(:,3) Maps real particle index to the ordered set
    
    tau = 0.2*(2*R/cells)/v0; % Timestep
    coeff = 0.5*eff_num*pi*diam*diam*tau/(L*L*L/(cells^3)); % Not quite sure where it comes from
    
    vrmax = zeros(cells,cells,cells) + 3*v0; % Initial vmax
    selxtra = zeros(cells,cells,cells); % The 'Rem' part in N_pairs on page 4 in the GPU article
    coltot = 0; % Total number of collisions
    oldColTot = 0;
    
    numParticles = 0;
    
    % Particle arrays
    p = zeros(100*particles,1);
    up = linspace(1,100*particles,100*particles);
    r = zeros(100*particles,3);
    v = zeros(100*particles,3);
    
    sprintf('Starting simulation with')
    sprintf('Particles: %d',particles)
    sprintf('Cells (each dimension): %d',cells)
    sprintf('Cells (total): %d',cells^3)
    sprintf('v0 = %f m/s',v0)
    sprintf('T = %f kelvin',T)
    sprintf('tau = %f ns',tau * 10^9)
    sprintf('Each particle represents %f atoms',eff_num)
    
    E = zeros(1,timesteps);
    
    totalNumberOfParticles = zeros(timesteps,1);
    
    % [numParticles,p,up,r,v] = createParticles(numParticles,p,up,r,v,R,sigma,500);
    
    E(1) = 0.5*sum(sum(v.^2,2));
    sprintf('Energy at t=0: %f J',E(1))
    sprintf('Mean velocity = %f m/s',mean(mean(abs(v),2)))
    
    animate = true;
    if(animate)
        fig1=figure(1);
        winsize = get(fig1,'Position');
        winsize(1:2) = [0 0];
        A=moviein(timesteps,fig1,winsize);
        set(fig1,'NextPlot','replacechildren');
    end
    
    t = linspace(0,2*pi,100);
    circleX = R*cos(t);
    circleY = R*sin(t);
    
    for i=1:timesteps
       totalNumberOfParticles(i) = numParticles;
        
       if(energyVsTime) E(i) = 0.5*sum(sum(v.^2,2)); end
       
       sortData = sortParticles(r,L,cells,sortData,numParticles,p); % Put particles in cells
       
       [r,v] = mover(r,cells,v,tau,L,R,sigma,numParticles,p); % Move and collide with walls
       [col, v, vrmax, selxtra] = collide(v,vrmax,selxtra,coeff,sortData,cells,L,numParticles,p); % Collide with other particles
       
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,100) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           
           sprintf('Done %d of %d steps. %d collisions (%d new)',i,timesteps,coltot,dCol)
           sprintf('Total energy: %f',E(i))
           % plotMeanXVelocity(r,v,R,numParticles,p)
       end
       
       
       % Clear figure and plot this timestep
        particleIndices = p(1:numParticles);
           figure(1);
           clf(fig1);
           plot(r(particleIndices,1),r(particleIndices,2),'d');
           hold on;
           % plot(circleX,circleY,'g');
           xlabel('z [m]');
           ylabel('y [m]');
           axis('equal')
           
           axis([0 L -R R]);

           title(sprintf('t = %f ns',10^9*i*tau));
           A(:,i)=getframe(fig1,winsize); % save to movie
           
           [r,v,p,up,numParticles] = cleanUpParticles(p,up,numParticles,r,v,L);
           [numParticles,p,up,r,v] = createParticles(numParticles,p,up,r,v,R,sigma,10);
    end
    
    if(energyVsTime)
        t = linspace(0,timesteps*tau,timesteps);
        plotEnergy(t,E);
        axis([0 max(t) 0 max(E)]);
        axis 'auto x';
        
        % figure
        % plot(t,totalNumberOfParticles);
    end
    
end