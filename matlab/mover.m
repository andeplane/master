function [r,v] = mover(r,sd,cells,v,tau,L,mpv)
    yold = r(:,2);
    r = r + v*tau;
    r(:,1) = mod(r(:,1)+1000*L,L); %Periodic boundary conditions
    [r,v] = collideWithWalls(r,sd,cells,v,tau,L,mpv,yold);
    
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

function [r,v] = collideWithWalls(r,sd,cells,v,tau,L,mpv,yold)
    cj = @(k) ceil(k/cells);
    
    direction = [1, -1];
    ywall = [0, L];
    stddev = mpv/sqrt(2);
              
    for jcell=1:cells*cells
       j = cj(jcell);
       
       if j == 1 || j == cells % close to the walls
           particlesInCell = sd(jcell,1);
           
           % Loop through all particles and see if they are colliding
           for ipart = 1:particlesInCell; 
               
              ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index
              flag = 0;
              
              if(r(ip1,2) <= 0) flag = 1; end
              if(r(ip1,2) >= L) flag = 2; end
              
              if(flag > 0)
                % we collided
                v(ip1,2) = direction(flag)*sqrt(-log(1-rand())) * mpv;
                v(ip1,1) = stddev*normrnd(0,1);

                dtr = tau*(r(ip1,2)-ywall(flag))/(r(ip1,2)-yold(ip1));
                
                r(ip1,2) = ywall(flag) + v(ip1,2)*dtr;
              end
           end
       end
    end
end