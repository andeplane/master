function particles = createParticles(particles,sigma,newParticlesCount,dt,xrange,yrange,zrange)
    newParticleIndices = zeros(1,newParticlesCount);
    
    for i=1:newParticlesCount
       particles.count = particles.count+1;
       particleIndex = particles.unused(particles.count);
       
       particles.indices(particles.count) = particleIndex;
       particles.rtop(particleIndex) = particles.count;
       
       newParticleIndices(i) = particleIndex;
       
       particles.r(particleIndex,:) = [xrange(1)+rand()*(xrange(2) - xrange(1)),yrange(1)+rand()*(yrange(2) - yrange(1)),zrange(1)+rand()*(zrange(2) - zrange(1))];
       
       particles.v(particleIndex,:) = normrnd(0,sigma,1,3); % Maxwell distribution
       particles.v(particleIndex,1) = particles.v(particleIndex,1) + 3*sigma;
    end
    
    [tetids, bcs] = pointLocation(dt, particles.r(newParticleIndices,:));
    sprintf('Number of particles outside before: %d',sum(isnan(tetids)))
    particlesOutside = find(isnan(tetids));
    particles = removeParticles(particles,particlesOutside,newParticleIndices);
    
    [tetids, bcs] = pointLocation(dt, particles.r(newParticleIndices,:));
    sprintf('Number of particles outside after: %d',sum(isnan(tetids)))
    
    % [r,v,p,up,numParticles] = removeParticles(p,up,numParticles,r,v,particlesOutside,newParticleIndices);
    % sprintf('Removing %d particles that are outside',length(particlesOutside))
end

% function [r,v,p,up,numParticles] = removeParticles(p,up,numParticles,r,v,indices,newParticleIndices)
function particles = removeParticles(particles,indices,newParticleIndices)
    for i=1:length(indices)
        particleIndex = newParticleIndices(indices(i)); % This is indices in the newParticleIndices-array
        pidx = particles.rtop(particleIndex);
        
        % sprintf('Working on particle %d/%d with r-index=%d and pidx=%d',i,length(indices),particleIndex,pidx)
        
        particles.unused(particles.count) = particleIndex; % This particle is now unused
        particles.indices(pidx) = 0; % Removed
        
        particles.count = particles.count-1;
        
        particles.v(particleIndex,:) = [0 0 0];
        particles.r(particleIndex,:) = [0 0 0];
    end
    
    newIndices = find(particles.indices > 0);
    particles.indices(1:length(newIndices)) = newIndices;
end