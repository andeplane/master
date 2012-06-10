function [r,v] = mover(r,sd,cells,v,tau,L,R,stddev,numParticles,p)
    particleIndices = p(1:numParticles);
    rold = r;
    r(particleIndices,:) = r(particleIndices,:) + v(particleIndices,:)*tau;
    
    r(particleIndices,1) = mod(r(particleIndices,1)+1000*L,L); %Periodic boundary conditions
    
    [r,v] = collideWithWalls(r,cells,v,R,stddev,rold,p,numParticles);
end

function [r,v] = collideWithWalls(r,cells,v,R,stddev,rold,p,numParticles)
    cj = @(l) ceil((mod(l-1,cells*cells)+1)/cells);
    an = @(x,y) mod((atan2(y,x) + 100*pi),2*pi); % Calculate angle
    transform = @(vector,angle) ([cos(angle) -sin(angle); sin(angle) cos(angle)]*vector')';
    
    factor = sqrt(2)*stddev;
    
    for ipart=1:numParticles
       particleIndex = p(ipart); %Real particle index
       if(norm(r(particleIndex,2:3)) > R)
          % Collided with wall
          angle = an(r(particleIndex,3),r(particleIndex,2));
          
          r(particleIndex,:) = rold(particleIndex,:);
          
          % [z,y] are rotated
          vz = -sqrt(-log(rand())) * factor;
          vy = stddev*normrnd(0,1);
          vx = stddev*normrnd(0,1);
          
          vzy = transform([vz,vy],angle); % Rotate back to normal coordinates
          
          v(particleIndex,:) = [vx vzy([2,1])];
       end
       
       if(norm(r(particleIndex,2:3)) > R)
          sprintf('We have a particle outside :/') 
       end
    end
              
%     for jcell=1:cells*cells
%        j = cj(jcell);
%        
%        if j == 1 || j == cells % close to the walls
%            particlesInCell = sd(jcell,1);
%            
%            % Loop through all particles and see if they are colliding
%            for ipart = 1:particlesInCell; 
%               ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index in particle index list, jeez
%               particleIndex = p(ip1);
%               
%               flag = 0;
%               
%               if(r(particleIndex,2) <= 0) flag = 1; end
%               if(r(particleIndex,2) >= L) flag = 2; end
%               
%               if(flag > 0)
%                 % we collided
%                 
%                 v(particleIndex,2) = direction(flag)*sqrt(-log(1-rand())) * factor;
%                 v(particleIndex,1) = stddev*normrnd(0,1);
%                 
%                 dtr = tau*(r(particleIndex,2)-ywall(flag))/(r(particleIndex,2)-yold(particleIndex));
%                 
%                 r(particleIndex,2) = ywall(flag) + v(particleIndex,2)*dtr;
%               end
%            end
%        end
%     end
end