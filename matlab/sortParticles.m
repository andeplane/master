function sd = sortParticles(r,L,cells,sd,numParticles,p)
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells)+1;
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*cells + i;
    
    jx = zeros(numParticles,1);
    jy = zeros(numParticles,1);
    
    sd(:,1) = 0; % Reset the number of particles in each cell
    
    for ipart=1:numParticles
        %Place particle i in the correct cell based on its position. 
        particleIndex = p(ipart);
        
        i = ceil(r(particleIndex,1)*cells/L);
        j = ceil(r(particleIndex,2)*cells/L);
        
        if(r(particleIndex,2) > 2*L || r(particleIndex,2) < -2*L)
           sprintf('We have a particle (%d) that is at %f L, j=%d',ipart,r(ipart,2)/L,j) 
        end
        
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
    for ipart=1:numParticles
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