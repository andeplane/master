function dsmc(particles,timesteps,particlesPerCell)
    format long;
    
    % Plots
    animate = true;
    meanXVelocity = true;
    energyVsTime = true;
    plotPressure = false;
    
    % Calculate the required number of cells
    cells = ceil(sqrt(particles/particlesPerCell));
    
    % Some constants
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    mass = 6.63e-26;       % Mass of argon atom (kg)
    massgrams = mass*1000;
    
    diam = 3.66e-10;       % Effective diameter of argon atom (m)
    T = 273;               % Temperature (K)
    density = 1.78;        % Density of argon at STP (kg/m^3)
    L = 1e-6;              % System size is one micron
    
    eff_num = density/mass*L*L*L/particles;
    v0 = sqrt(3*boltz*T/mass); % Used to calculate initial vmax
    sigma = sqrt(boltz*T/mass); % Standard deviation
    
    % For small systems (testing) the number of cells^2 might be larger
    % than the number of particles
    dim = max(particles,cells*cells) + 10;
    
    sortData = zeros(dim,3);
    
    %sortData explanation:
       %cell_n = sortData(:,1) Number of particles in a specific cell
       %index = sortData(:,2) Ordered set of particles
       %Xref = sortData(:,3) Maps real particle index to the ordered set
    
    tau = 0.2*(L/cells)/v0; % Timestep
    coeff = 0.5*eff_num*pi*diam*diam*tau/(L*L*L/(cells*cells)); % Not quite sure where it comes from
    vrmax = zeros(cells,cells) + 3*v0; % Initial vmax
    selxtra = zeros(cells,cells); % The 'Rem' part in N_pairs on page 4 in the GPU article
    coltot = 0; % Total number of collisions
    oldColTot = 0;
    
    % Prepare movie
    if(animate)
        fig1=figure(1);
        winsize = get(fig1,'Position');
        winsize(1:2) = [0 0];
        A=moviein(timesteps,fig1,winsize);
        set(fig1,'NextPlot','replacechildren');
    end
    
    xmax = L*1.05;
    xmin = L-xmax;
    
    r = L*rand(particles,3); % Initiate random positions
    v = normrnd(0,sigma,particles,3); % Maxwell distribution
    
    v(:,1) = v(:,1) + 3*sigma; % add gas velocity in x-direction
    
    % set z and v_z to 0
    v(:,3) = 0; 
    r(:,3) = 0;
    
    % Calculate pressure (ideal gas)
    % 39.948 g/mol
    chemAmount = particles*eff_num*massgrams / 39.948;
    pTot = chemAmount*208*T/(L*L); % total pressure
    
    sprintf('Starting simulation with')
    sprintf('Particles: %d',particles)
    sprintf('Cells (each dimension): %d',cells)
    sprintf('Cells (total): %d',cells*cells)
    sprintf('v0 = %f m/s',v0)
    sprintf('T = %f kelvin',T)
    sprintf('P = %f Pa',pTot)
    sprintf('tau = %f ns',tau * 10^9)
    sprintf('Each particle represents %f atoms',eff_num)
    
    E = zeros(1,timesteps);
    Pmax = 0;
    
    for i=1:timesteps
        if(energyVsTime) E(i) = 0.5*sum(sum(v.^2,2)); end
        
       if(animate)
           % Clear figure and plot this timestep
           figure(1);
           clf(fig1);
           plot(r(:,1),r(:,2),'d');
           hold on;
           plot([0 L],[L,L],'r');
           plot([0 L],[0,0],'r');
           plot([L L],[0,L],'c');
           plot([0 0],[0,L],'c');
           xlabel('x [m]');
           ylabel('y [m]');
           axis('equal')
           % axis([xmin xmax xmin xmax]);

           title(sprintf('t = %f ns',10^9*i*tau));
           A(:,i)=getframe(fig1,winsize); % save to movie
       end
       v_before = v;
       
       sortData = sortParticles(r,L,cells,particles,sortData); %Put particles in cells
       [r,v] = mover(r,sortData,cells,v,tau,L,v0); % Move and collide with walls
       [col, v, vrmax, selxtra] = collide(v,vrmax,selxtra,coeff,sortData,cells,L); % Collide with other particles
       deltaV = v - v_before;
       
       [P,Pmax] = calculatePressure(sortData,cells,v,mass,L,plotPressure,Pmax,eff_num,deltaV);
       
       
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,100) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           if meanXVelocity plotMeanXVelocity(r,v,L,particles); end
           
           
           sprintf('Done %d of %d steps. %d collisions (%d new)',i,timesteps,coltot,dCol)
       end
    end
    
    if(energyVsTime)
        t = linspace(0,timesteps*tau,timesteps);
        plotEnergy(t,E);
    end
    
end