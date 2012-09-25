function dsmc(timesteps)
    close all;
    format long;
    filename = 'shrinktube.obj';
    
    [dt, points, fn] = triangulateObjFile(filename);
    xrange = [min(points(:,1)) max(points(:,1))]
    yrange = [min(points(:,2)) max(points(:,2))]
    zrange = [min(points(:,3)) max(points(:,3))]
    
    particlesPerCell = 25;
    particleCount = 20000;
    
    % Calculate the required number of cells
    cells = ceil((particleCount/particlesPerCell)^(1/3));
    
    % Some constants
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    mass = 6.63e-26;       % Mass of argon atom (kg)
    diam = 3.66e-10;       % Effective diameter of argon atom (m)
    T = 273;               % Temperature (K)
    density = 1.78;        % Density of argon at STP (kg/m^3)
    R = 100*1e-9;              % Tube radius
    L = 1000*1e-9;               % Tube length
    
    volume = pi*R*R*L;
    
    eff_num = density/mass*volume/particleCount;
    v0 = sqrt(3*boltz*T/mass); % Used to calculate initial vmax
    sigma = sqrt(boltz*T/mass); % Standard deviation
    
    % For small systems (testing) the number of cells^2 might be larger
    % than the number of particles
    dim = max(100*particleCount,cells^3) + 10;
    sortData = zeros(dim,3);
    
    %sortData explanation:
       %cell_n = sortData(:,1) Number of particles in a specific cell
       %index = sortData(:,2) Ordered set of particles
       %Xref = sortData(:,3) Maps real particle index to the ordered set
    
    tau = 0.05  *(2*R/cells)/v0; % Timestep
    
    vrmax = zeros(cells,cells,cells) + 3*v0; % Initial vmax
    selxtra = zeros(cells,cells,cells); % The 'Rem' part in N_pairs on page 4 in the GPU article
    coltot = 0; % Total number of collisions
    oldColTot = 0;
    
    numParticles = 0;
    
    % Particle arrays
    p = zeros(100*particleCount,1);
    rtop = zeros(100*particleCount,1);
    up = linspace(1,100*particleCount,100*particleCount);
    r = zeros(100*particleCount,3);
    v = zeros(100*particleCount,3);
    
    sprintf('Starting simulation with')
    sprintf('Particles: %d',particleCount)
    sprintf('Cells (each dimension): %d',cells)
    sprintf('Cells (total): %d',cells^3)
    sprintf('v0 = %f m/s',v0)
    sprintf('T = %f kelvin',T)
    sprintf('tau = %f ns',tau * 10^9)
    sprintf('Each particle represents %f atoms',eff_num)
    
    totalNumberOfParticles = zeros(timesteps,1);
    
    particles = struct('count',numParticles,'indices',p,'rtop',rtop,'unused',up,'r',r,'v',v);
    
    particles = createParticles(particles,sigma,1/4*5000,dt,xrange,yrange,zrange);
    
    particleIndices = particles.indices(1:particles.count);
    
    sprintf('Mean velocity = %f m/s',mean(mean(abs(particles.v(particleIndices,:)),2)))
    
    fid = fopen('pos.xyz','w');
    
    for i=1:timesteps
       eff_num = density/mass*volume/particleCount;
       coeff = 0.5*eff_num*pi*diam*diam*tau/(volume/cells^3); % Not quite sure where it comes from
       totalNumberOfParticles(i) = particles.count;
       sortData = sortParticles(L,R,cells,sortData,particles); % Put particles in cells
       
       particles = mover(tau,dt,particles); % Move and collide with walls
       
       [col, particles, vrmax, selxtra] = collide(vrmax,selxtra,coeff,sortData,cells,particles); % Collide with other particles
       
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,10) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           
           sprintf('Done %d of %d steps. %d collisions (%d new)',i,timesteps,coltot,dCol)
           sprintf('Living particles: %f',particles.count)
       end
       % sprintf('Cleaned up particles')
       % [r,v,p,up,numParticles] = cleanUpParticles(p,up,numParticles,r,v,L,R);
       % newParticles = particles - numParticles;
       % [numParticles,p,up,r,v] = createParticles(numParticles,p,up,r,v,sigma,20,newParticles,dt,xrange,yrange,zrange);
       
       particleIndices = particles.indices(1:particles.count);
       fprintf(fid,'%d\n a\n',particles.count);
       fprintf(fid,'H %f %f %f\n',particles.r(particleIndices,:)'*10^6);
    end
    
    fclose(fid);
end