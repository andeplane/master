function [numParticles,p,up,r,v] = createParticles(numParticles,p,up,r,v,R,sigma,newParticlesCount)
    for i=1:newParticlesCount
       numParticles = numParticles + 1;
       particleIndex = up(numParticles);
       p(numParticles) = particleIndex;
       
       theta = 2*pi*rand();
       radial = R*sqrt(rand());
       z = radial*cos(theta);
       y = radial*sin(theta);
       
       r(particleIndex,:) = [0,y,z];
       v(particleIndex,:) = normrnd(0,sigma,1,3); % Maxwell distribution
       v(particleIndex,1) = v(particleIndex,1) + 3*sigma;
    end
end