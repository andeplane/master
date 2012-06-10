function [r,v] = mover(r,sd,cells,v,tau,L,stddev,numParticles,p)
    particleIndices = p(1:numParticles);
    yold = r(:,2);
    r(particleIndices,:) = r(particleIndices,:) + v(particleIndices,:)*tau;
    
    r(particleIndices,1) = mod(r(particleIndices,1)+1000*L,L); %Periodic boundary conditions
    
    [r,v] = collideWithWalls(r,sd,cells,v,tau,L,stddev,yold,p);
    
    %newLivingParticles = 0;
    %newParticleIndices = zeros(10^7,1);
    %
    %for ipart=1:livingParticles
    %    pid = particleIndices(ipart);
    %    if(r(pid,1) <= L && r(pid,1) >= 0) 
    %        newLivingParticles = newLivingParticles + 1;
    %        newParticleIndices(newLivingParticles) = pid;
    %    end
    %end
end

function [r,v] = collideWithWalls(r,sd,cells,v,tau,L,stddev,yold,p)
    cj = @(k) ceil(k/cells);
    boltz = 1.3806e-23;    % Boltzmann (J/K)
    
    direction = [1, -1];
    ywall = [0, L];
    factor = sqrt(2)*stddev;
              
    for jcell=1:cells*cells
       j = cj(jcell);
       
       if j == 1 || j == cells % close to the walls
           particlesInCell = sd(jcell,1);
           
           % Loop through all particles and see if they are colliding
           for ipart = 1:particlesInCell; 
              ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index in particle index list, jeez
              particleIndex = p(ip1);
              
              flag = 0;
              
              if(r(particleIndex,2) <= 0) flag = 1; end
              if(r(particleIndex,2) >= L) flag = 2; end
              
              if(flag > 0)
                % we collided
                
                v(particleIndex,2) = direction(flag)*sqrt(-log(1-rand())) * factor;
                v(particleIndex,1) = stddev*normrnd(0,1);
                
                dtr = tau*(r(particleIndex,2)-ywall(flag))/(r(particleIndex,2)-yold(particleIndex));
                
                r(particleIndex,2) = ywall(flag) + v(particleIndex,2)*dtr;
              end
           end
       end
    end
end