function sd = sortParticles(L,R,cells,sd,particles)
% function sd = sortParticles(r,L,R,cells,sd,numParticles,p)
% sortData = sortParticles(L,R,cells,sortData,particles);
    %cell number as function of i, j
    cl = @(i,j,k) (k-1)*cells^2 + (j-1)*cells + i;
    
    jx = zeros(particles.count,1);
    jy = zeros(particles.count,1);
    jz = zeros(particles.count,1);
    
    sd(:,1) = 0; % Reset the number of particles in each cell
    
    for ipart=1:particles.count
        %Place particle i in the correct cell based on its position. 
        particleIndex = particles.indices(ipart);
        
        rx = particles.r(particleIndex,1);
        ry = particles.r(particleIndex,2);
        rz = particles.r(particleIndex,3);
        
        i = ceil(rx*cells/L);
        j = ceil((ry+R)*cells/(2*R));
        k = ceil((rz+R)*cells/(2*R));
        
        if(i < 1) i = 1; end
        if(i > cells) i = cells; end
        if(j < 1) j = 1; end
        if(j > cells) j = cells; end
        if(k > cells) k = cells; end
        if(k < 1) k = 1; end
        
        jx(ipart) = i;
        jy(ipart) = j;
        jz(ipart) = k;
        
        % This is cell number k, calculate with above eqs
        l = cl(i,j,k);
        
        sd(l,1) = sd(l,1) + 1;
    end
    
    % Create an ordered set of particles, sorted by cell index
    m = 1;
    for l=1:cells^3
       sd(l,2) = m;
       m = m + sd(l,1);
    end
    
    % Create a mapping between the sorted set to the real particle indices
    temp = zeros(cells^3,1);
    for ipart=1:particles.count
        i = jx(ipart);
        j = jy(ipart);
        k = jz(ipart);
        
        % cell k
        l = cl(i,j,k);
        if(l>cells*cells*cells)
           sprintf('l (%d) too large, i,j=%d,%d',k,i,j) 
        end
        
        m = sd(l,2) + temp(l);
        sd(m,3) = ipart;
        temp(l) = temp(l) + 1;
    end
end