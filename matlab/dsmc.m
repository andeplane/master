function dsmc(particles,timesteps,cells,animate)
    format long;
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells);
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*4 + i;
    
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
       % [r,v] = collideWithCircle(r,sortData,cells,v,tau,radius,xc,yc);
       
       
       r(:,1) = mod(r(:,1)+1000*L0,L0); %Periodic boundary conditions
       % r = mod(r+1000*L,L); %Periodic boundary conditions
       
       if(animate)
           % Clear figure and plot this timestep
           clf(fig1);
           plot(r(:,1),r(:,2),'d');
           hold on;
           plot([0 L0],[L,L],'r');
           plot([0 L0],[0,0],'r');
           plot([L0 L0],[0,L],'c');
           plot([0 0],[0,L],'c');
           plot(circleX,circleY,'g');
           axis('equal')
           axis([xmin xmax xmin xmax]);

           title(sprintf('t = %f ns',1000000*i*tau));
           A(:,i)=getframe(fig1,winsize); % save to movie
       end
       
       sortData = sort(r,L,cells,particles,sortData); %Put particles in cells
       
       [r,v] = collideWithEnv(r,sortData,cells,v,tau,L);
       
       
       [col, v, vrmax, selxtra] = collide(v,vrmax,selxtra,coeff,sortData,cells);
       
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,100) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           
           %k = 1000;
           
           %xvel = zeros(k,1);
           %count = zeros(k,1);
           %for ipart=1:particles
           %    verCell = ceil(k*r(ipart,2)/L);
           %    if(verCell > k) verCell = k; end
           %    xvel(verCell) = xvel(verCell) + v(ipart,1);
           %    count(verCell) = count(verCell)+1;
           %end
           
           %figure(2)
           %xvel = smooth(xvel./count);
           %plot(xvel)
           
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

function [r,v] = collideWithEnv(r,sd,cells,v,tau,L)
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells);
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*4 + i;
    friction = 0.1;
    for jcell=1:cells*cells
       j = cj(jcell);
       
       %if j == 1 || j == cells % close to the walls
           particlesInCell = sd(jcell,1);
           
           % Loop through all particles and see if they are colliding
           for ipart = 1:particlesInCell; 
               
              ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index
              
              if r(ip1,2) < 0
                 oldR = r(ip1,:) - v(ip1,:)*tau; %Go back in time
                 dt = oldR(2) / v(ip1,2); % Calculate time until collision
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*dt; % Move the particles less
                 v(ip1,2) = -v(ip1,2); % Invert velocity
                 v(ip1,1) = v(ip1,1)*friction;
                 % Move the particles the rest of the time
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*(tau - dt);
              end
              
              if r(ip1,2) > L
                 oldR = r(ip1,:) - v(ip1,:)*tau; %Go back in time
                 dt = (L - oldR(2)) / v(ip1,2); % Calculate time until collision
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*dt; % Move the particles less
                 v(ip1,2) = -v(ip1,2); % Invert velocity
                 v(ip1,1) = v(ip1,1)*friction;
                 % Move the particles the rest of the time
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*(tau - dt);
              end
              
           end
       %end
    end
end

function [r,v] = collideWithCircle(r,sd,cells,v,tau,radius,xc,yc)
    for jcell=1:cells*cells
       %if j == 1 || j == cells % close to the walls
           particlesInCell = sd(jcell,1);
            
           % Loop through all particles and see if they are colliding
           for ipart = 1:particlesInCell; 
               
              ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index
              
              if(norm(r(ip1,:) - [xc,yc,0]) < radius)
                 %collided
                 
                 oldR = r(ip1,:) - v(ip1,:)*tau; %Go back in time
                 % Calculate time of colision
                 x0 = oldR(1);
                 y0 = oldR(2);
                 vx = v(ip1,1);
                 vy = v(ip1,2);
                 
                 a = v(ip1,1)*v(ip1,1) + v(ip1,2)*v(ip1,2);
                 b = 2*(x0*vx + y0*vy - xc*vx - yc*vy);
                 c = -radius^2 - 2*xc*x0 - 2*yc*y0 + x0^2 + y0^2 + yc^2 + xc^2;
                 t0 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
                 t1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
                 
                 t = min(t0,t1);
                 
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*t;
                 v(ip1,:) = -v(ip1,:);
                 
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*(tau-t);
              end
              
           end
       %end
    end

end

function [col, v, crmax, selxtra] = collide(v,crmax,selxtra,coeff,sd,cells) 
    col = 0;
    for jcell = 1:cells*cells
        number = sd(jcell,1); %Number of particles in this cell
        if number < 2
            continue; 
        end
       
        % How many collisions (N_pairs in gpu article)
        select = coeff*number*(number-1)*crmax(jcell) + selxtra(jcell);
        nsel = floor(select);
        selxtra(jcell) = select - nsel; % Rem
        crm = crmax(jcell); %Max relative speed in this cell
        
        for isel=1:nsel
           t1 = floor(rand()*number); % Random index in this cell
           t2 = mod((floor(t1+rand()*(number-1))+1),number); % Another random index
           
           ip1 = sd(t1+sd(jcell,2),3); % Actual particle index
           ip2 = sd(t2+sd(jcell,2),3);
           
           % Calculate relative speed
           cr = sqrt((v(ip1,1) - v(ip2,1))^2 + (v(ip1,2) - v(ip2,2))^2 + (v(ip1,3) - v(ip2,3))^2);
           if cr > crm
               % Update max relative speed
               crm = cr;
           end
           
           if cr/crmax(jcell) > rand()
              % We have a collision
              col = col + 1; 
              
              % Center of mass velocity
              vcm = 0.5*(v(ip1,:) + v(ip2,:));
              
              cos_th = 1 - 2*rand(); % Random cosine
              sin_th = sqrt(1 - cos_th^2);
              phi = 2*pi*rand();
              vrel = zeros(1,3);
              vrel(1) = cr*cos_th;
              vrel(2) = cr*sin_th*cos(phi);
              vrel(3) = cr*sin_th*sin(phi);
              
              v(ip1,:) = vcm(:) + 0.5*vrel(:);
              v(ip2,:) = vcm(:) - 0.5*vrel(:);
           end
           
           crmax(jcell) = crm;
        end
    end
end

function sd = sort(r,L,cells,particles,sd)
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells);
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*4 + i;
    
    jx = zeros(particles,1);
    jy = zeros(particles,1);
    
    sd(:,1) = 0; % Reset the number of particles in each cell
    
    for ipart=1:particles
        %Place particle i in the correct cell based on its position. 
        i = ceil(r(ipart,1)*cells/L);
        j = ceil(r(ipart,2)*cells/L);
        
        if(i < 1) i = 1; end
        if(i > cells) i = cells; end
        if(j < 1) j = 1; end
        if(j > cells) j = cells; end
        
        jx(ipart) = i;
        jy(ipart) = j;
        
        % This is cell number k, calculate with above eqs
        k = ck(i,j);
        
        sd(k,1) = sd(k,1) + 1;
    end
    
    % Create an ordered set of particles, sorted by cell index
    m = 1;
    for k=1:cells*cells
       sd(k,2) = m;
       m = m + sd(k,1);
    end
    
    % Create a mapping between the sorted set to the real particle indices
    temp = zeros(cells*cells,1);
    for ipart=1:particles
        i = jx(ipart);
        j = jy(ipart);
        
        % cell k
        k = ck(i,j);
        if(k>cells*cells)
           sprintf('k (%d) too large, i,j=%d,%d',k,i,j) 
        end
        
        l = sd(k,2) + temp(k);
        sd(l,3) = ipart;
        temp(k) = temp(k) + 1;
    end
end