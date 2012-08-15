function dsmc(timesteps)
    close all;
    format long;
    
    % Plots
    doAnimate = true;
    energyVsTime = false;
    
    particlesPerCell = 25;
    particles = 20000;
    
    % Calculate the required number of cells
    cells = ceil((particles/particlesPerCell)^(1/3));
    
    % Some constants
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    mass = 6.63e-26;       % Mass of argon atom (kg)
    diam = 3.66e-10;       % Effective diameter of argon atom (m)
    T = 273;               % Temperature (K)
    density = 1.78;        % Density of argon at STP (kg/m^3)
    R = 5e-7;              % Tube radius
    L = 5*R;              % Tube length
    
    
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
    
    [numParticles,p,up,r,v] = createParticles(numParticles,p,up,r,v,R,L,sigma,1/4*5000,true);
    
    E(1) = 0.5*sum(sum(v.^2,2));
    sprintf('Energy at t=0: %f J',E(1))
    particleIndices = p(1:numParticles);
    
    sprintf('Mean velocity = %f m/s',mean(mean(abs(v(particleIndices,:)),2)))
    
    if(doAnimate)
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
       eff_num = density/mass*volume/particles;
       coeff = 0.5*eff_num*pi*diam*diam*tau/(volume/cells^3); % Not quite sure where it comes from
       totalNumberOfParticles(i) = numParticles;
        
       if(energyVsTime) E(i) = 0.5*sum(sum(v.^2,2)); end
       
       sortData = sortParticles(r,L,R,cells,sortData,numParticles,p); % Put particles in cells
       % sprintf('Moving...')
       [r,v] = mover(r,cells,v,tau,R,L,sigma,numParticles,p); % Move and collide with walls
       % sprintf('Moved, done. Colliding...')
       [col, v, vrmax, selxtra] = collide(v,vrmax,selxtra,coeff,sortData,cells,p); % Collide with other particles
       % sprintf('Collided, done. Cleaning up...')
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,10) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           
           % plotMeanXVelocity(r,v,R,numParticles,p)
           
           sprintf('Done %d of %d steps. %d collisions (%d new)',i,timesteps,coltot,dCol)
           sprintf('Total energy: %f',E(i))
           sprintf('Living particles: %f',numParticles)
       end
       % sprintf('Cleaned up particles')
       [r,v,p,up,numParticles] = cleanUpParticles(p,up,numParticles,r,v,L,R);
       newParticles = particles - numParticles;
       [numParticles,p,up,r,v] = createParticles(numParticles,p,up,r,v,R,L,sigma,20,newParticles);
       
       if(doAnimate) animate(numParticles,p,r,R,L,i,tau,fig1,circleX,circleY,winsize,A); end
    end
    
    if(energyVsTime)
        t = linspace(0,timesteps*tau*10^9,timesteps);
        plotEnergy(t,E);
        axis([0 max(t) 0 max(E)]);
        axis 'auto x';
        xlabel('Time [ns]')
        xlabel('Energy')
        
        % figure
        % plot(t,totalNumberOfParticles);
    end
    
end