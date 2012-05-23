function dsmc(particles,timesteps,cells)
    format long;
    
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    mass = 6.63e-26;       % Mass of argon atom (kg)
    diam = 3.66e-10;       % Effective diameter of argon atom (m)
    T = 273;               % Temperature (K)
    density = 1.78;        % Density of argon at STP (kg/m^3)
    L = 1e-6;              % System size is one micron
    
    eff_num = density/mass*L*L*L/particles;
    sprintf('Each particle represents %f atoms',eff_num)
    v0 = sqrt(3*boltz*T/mass);
    
    x = zeros(particles);
    v = zeros(particles,3);
    sigma = sqrt(boltz*T/mass);
    mu = sqrt(2*boltz*T/mass);
    
    sortData = zeros(particles,3);
    
    %cell_n = sortData(:,1) cells long
    %index = sortData(:,2) cells long
    %Xref = sortData(:,3) particles long
    
    for i=1:particles % does matlab include particles?
        x(i) = L*rand();
        plusMinus = 2*round(rand()) - 1;
        %v(i,1) = normrnd(0,sigma);
        v(i,1) = v0*plusMinus;
    end
    
    vmagI = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    
    tau = 0.2*(L/cells)/v0
    vrmax = zeros(cells) + 3*v0;
    selxtra = zeros(cells);
    coeff = 0.5*eff_num*pi*diam*diam*tau/(L*L*L/cells); % Hva er dette?
    coltot = 0;
    
    for i=1:timesteps
       for n=1:particles
          x(n) = x(n) + v(n,1)*tau;
          x(n) = mod(x(n)+L,L); % Periodic boundaries
       end
       
       sortData = sort(x,L,cells,particles,sortData);
       %int col = colider(v,vrmax,tau,seed,selxtra,coeff,sortData);
       %coltot += col;	// Increment collision count
       [col, v, vrmax, selxtra] = collide(v,vrmax,tau,selxtra,coeff,sortData,cells,particles);
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,100) < 1
          sprintf('Done %d of %d steps. %d collisions',i,timesteps,coltot)
       end
    end
    
    vmagF = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    
    vbin = 50:100:1050;    % Bins for histogram
    hist(vmagI,vbin);  title('Initial speed distribution');
    xlabel('Speed (m/s)');  ylabel('Number');
    
    %* Plot the histogram of the final speed distribution
    figure(2); clf;
    hist(vmagF,vbin);  
    title(sprintf('Final speed distribution'));
    xlabel('Speed (m/s)');  ylabel('Number');
    

end

function [col, v, crmax, selxtra] = collide(v,crmax,tau,selxtra,coeff,sd,cells,particles) 
    col = 0;
    for jcell = 1:cells
        number = sd(jcell,1);
        if number < 2
            continue; 
        end
       
        select = coeff*number*(number-1)*crmax(jcell) + selxtra(jcell);
        nsel = floor(select);
        selxtra(jcell) = select - nsel;
        crm = crmax(jcell);
        for isel=1:nsel
           k = floor(rand()*number);
           kk = mod((floor(k+rand()*(number-1))+1),number); % number;
           ip1 = sd(k+sd(jcell,2),3);
           ip2 = sd(kk+sd(jcell,2),3);
           % Calculate relative speed
           cr = sqrt((v(ip1,1) - v(ip2,1))^2 + (v(ip1,2) - v(ip2,2))^2 + (v(ip1,3) - v(ip2,3))^2);
           if cr > crm
               % Update max relative speed
               crm = cr;
           end
           
           if cr/crmax(jcell) > rand()
               % we did collide
              col = col + 1; 
              vcm = zeros(1,3);
              for k=1:3
                 vcm(k) = 0.5*(v(ip1,k) + v(ip2,k));
              end
              cos_th = 1 - 2*rand();
              sin_th = sqrt(1 - cos_th^2);
              phi = 2*pi*rand();
              vrel = zeros(1,3);
              vrel(1) = cr*cos_th;
              vrel(2) = cr*sin_th*cos(phi);
              vrel(3) = cr*sin_th*sin(phi);
              for k=1:3
                 v(ip1,k) = vcm(k) + 0.5*vrel(k);
                 v(ip2,k) = vcm(k) - 0.5*vrel(k);
              end
           end
           
           crmax(jcell) = crm;
        end
    end
end

function sd = sort(x,L,cells,particles,sd)
    jx = zeros(particles,1);
    
    sd(:,1) = 0; % Reset the number of particles in each cell
    
    for i=1:particles
        %Place particle i in the correct cell based on its position. 
        
        j = ceil(x(i)*cells/L);
        jx(i) = j;
        sd(j,1) = sd(j,1) + 1;
    end
    
    % Create a set order so that the particles are ordered so the
    % first N particles are the N particles in cell 1, the next M particles
    % are the M particles in cell 2 etc
    m = 1;
    for i=1:cells
       sd(i,2) = m;
       m = m + sd(i,1);
    end
    
    % Create a mapping between the sorted set to the real particle indices
    temp = zeros(cells,1);
    for i=1:particles
       jcell = jx(i);
       k = sd(jcell,2) + temp(jcell);
       sd(k,3) = i;
       temp(jcell) = temp(jcell) + 1;
    end
end