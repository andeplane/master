function [col, v, crmax, selxtra] = collide(v,crmax,selxtra,coeff,sd,cells,p) 
    ci = @(l) mod(l-1,cells)+1;
    cj = @(l) ceil((mod(l-1,cells*cells)+1)/cells);
    ck = @(l) ceil(l/cells^2);
    
    col = 0;
    colCells = zeros(cells,cells);
    
    for jcell = 1:cells^3
        number = sd(jcell,1); %Number of particles in this cell
        
        % sprintf('Colliding cell %d, partciles: %d',jcell,number)
        
        if number < 2
            continue; 
        end
        i = ci(jcell);
        j = cj(jcell);
        k = ck(jcell);
        
        % How many collisions (N_pairs in gpu article)
        
        select = coeff*number*(number-1)*crmax(jcell) + selxtra(jcell);
        
        nsel = floor(select);
        
        % sprintf('Colliding cell %d with %d particles and %d candidates',jcell,nsel, number)
        
        selxtra(jcell) = select - nsel; % Rem
        crm = crmax(jcell); %Max relative speed in this cell
        
        for isel=1:nsel
           t1 = floor(rand()*number); % Random index in this cell
           t2 = mod((floor(t1+rand()*(number-1))+1),number); % Another random index
           
           ip1 = sd(t1+sd(jcell,2),3); % Actual particle index
           ip2 = sd(t2+sd(jcell,2),3);
           
           particleIndex1 = p(ip1);
           particleIndex2 = p(ip2);
           
           % Calculate relative speed
           cr = norm(v(particleIndex1,:)-v(particleIndex2,:));
           
           % cr = sqrt((v(ip1,1) - v(ip2,1))^2 + (v(ip1,2) - v(ip2,2))^2 + (v(ip1,3) - v(ip2,3))^2);
           if cr > crm
               % Update max relative speed
               crm = cr;
           end
           
           if cr/crmax(jcell) > rand()
              % We have a collision
              col = col + 1;
              
              colCells(i,j) = colCells(i,j) + 1;
              
              % Center of mass velocity
              vcm = 0.5*(v(particleIndex1,:) + v(particleIndex2,:));
              
              cos_th = 1 - 2*rand(); % Random cosine
              sin_th = sqrt(1 - cos_th^2);
              phi = 2*pi*rand();
              vrel = zeros(1,3);
              vrel(1) = cr*cos_th;
              vrel(2) = cr*sin_th*cos(phi);
              vrel(3) = cr*sin_th*sin(phi);
              
              v(particleIndex1,:) = vcm(:) + 0.5*vrel(:);
              v(particleIndex2,:) = vcm(:) - 0.5*vrel(:);
           end
           
           crmax(jcell) = crm;
        end
    end
end