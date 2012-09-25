function [r,v,p,up,numParticles] = cleanUpParticles(p,up,numParticles,r,v,L,R)
    newP = zeros(length(p),1);
    newNumParticles = 0;
    oldNumParticles = numParticles;
    
    for ipart=1:oldNumParticles
       particleIndex = p(ipart);
       remove = false;
       % if(norm(r(particleIndex,2:3)) > R) remove = true; end
       if(r(particleIndex,1) > L) remove = true; end
       if(r(particleIndex,1) < 0) remove = true; end
       
       if(remove)
          up(numParticles) = particleIndex; % This particle is now unused
          numParticles = numParticles - 1;
          
          v(particleIndex,:) = [0 0 0];
          r(particleIndex,:) = [0 0 0];
       else
          newNumParticles = newNumParticles + 1;
          newP(newNumParticles) = particleIndex;
       end
    end
    
    p = newP;
    numParticles = newNumParticles;
end