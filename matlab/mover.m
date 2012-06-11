function [r,v] = mover(r,cells,v,tau,L,R,stddev,numParticles,p)
    particleIndices = p(1:numParticles);
    rold = r;
    r(particleIndices,:) = r(particleIndices,:) + v(particleIndices,:)*tau;
    
    % r(particleIndices,1) = mod(r(particleIndices,1)+1000*L,L); %Periodic boundary conditions
    
    [r,v] = collideWithWalls(r,cells,v,R,stddev,rold,p,numParticles);
end

function [r,v] = collideWithWalls(r,cells,v,R,stddev,rold,p,numParticles)
    cj = @(l) ceil((mod(l-1,cells*cells)+1)/cells);
    an = @(x,y) mod((atan2(y,x) + 4*pi),2*pi); % Calculate angle
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
end