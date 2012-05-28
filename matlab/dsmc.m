function dsmc(particles,timesteps,cells,animate)
    format long;
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells);
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*cells + i;
    
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    mass = 6.63e-26;       % Mass of argon atom (kg)
    diam = 3.66e-10;       % Effective diameter of argon atom (m)
    T = 273;               % Temperature (K)
    density = 1.78;        % Density of argon at STP (kg/m^3)
    L = 1e-6;              % System size is one micron
    
    eff_num = density/mass*L*L*L/particles;
    sprintf('Each particle represents %f atoms',eff_num)
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
    vrmax = zeros(cells,cells) + 3*v0; % Initial vmax
    selxtra = zeros(cells,cells); % The 'Rem' part in N_pairs on page 4 in the GPU article
    coltot = 0; % Total number of collisions
    oldColTot = 0;
    
    % Prepare movie
    fig1=figure(1);
    winsize = get(fig1,'Position');
    winsize(1:2) = [0 0];
    A=moviein(timesteps,fig1,winsize);
    set(fig1,'NextPlot','replacechildren');
    
    L0 = L;
    
    xmax = L0*1.05;
    xmin = L0-xmax;
    
    radius = L0/4;
    xc = L0/2;
    yc = L0/2;
    [circleX, circleY] = getCircle(xc,yc,radius);
    
    r = L*rand(particles,3); % Initiate random positions
    r(:,3) = 0;
    
    for ipart=1:particles
       dr = [xc,yc,0] - r(ipart,:);
       drLen = norm(dr);
       while(drLen < radius)
           r(ipart,1:2) = L*rand(1,2);
           
           dr = [xc,yc,0] - r(ipart,:);
           drLen = norm(dr);
       end
    end
    
    v = normrnd(0,sigma,particles,3); % Maxwell distribution
    v(:,1) = v(:,1) + 3*sigma;
    v(:,3) = 0; % set z to 0
    
    E0 = 0.5*v(:,:).^2;
    E = E0;
    
    for i=1:timesteps
       %L = max(L - L0*0.001,0.6*L0);
       
       coeff = 0.5*eff_num*pi*diam*diam*tau/(L*L*L/(cells*cells)); % Not quite sure where it comes from
        
       % Integrate position
       r = r + v*tau;
       
       r(:,1) = mod(r(:,1)+1000*L0,L0); %Periodic boundary conditions
       
       % r = mod(r+1000*L,L); %Periodic boundary conditions
       
       % [r,v] = collideWithCircle(r,sortData,cells,v,tau,radius,xc,yc,L,particles);
       
       if(animate)
           % Clear figure and plot this timestep
           figure(1);
           clf(fig1);
           plot(r(:,1),r(:,2),'d');
           hold on;
           plot([0 L0],[L,L],'r');
           plot([0 L0],[0,0],'r');
           plot([L0 L0],[0,L],'c');
           plot([0 0],[0,L],'c');
           plot(circleX,circleY,'g');
           axis('equal')
           % axis([xmin xmax xmin xmax]);

           title(sprintf('t = %f ns',1000000*i*tau));
           A(:,i)=getframe(fig1,winsize); % save to movie
       end
       
       sortData = sortParticles(r,L,cells,particles,sortData); %Put particles in cells
       
       [r,v] = collideWithEnv(r,sortData,cells,v,tau,L);
       
       [col, v, vrmax, selxtra] = collide(v,vrmax,selxtra,coeff,sortData,cells);
       
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,100) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           % plotMeanXVelocity(r,v,L,particles);
           
           sprintf('Done %d of %d steps. %d collisions (%d new)',i,timesteps,coltot,dCol)
       end
    end
    
    %Play movie
    movie(fig1,A,1,3,winsize);
end

function [x,y] = getCircle(xc,yc,radius)
    t = linspace(0,2*pi,30);
    x = radius*cos(t) + xc;
    y = radius*sin(t) + yc;
end