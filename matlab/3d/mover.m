function [r,v] = mover(r,cells,v,tau,R,L,stddev,numParticles,p)
    particleIndices = p(1:numParticles);
    rold = r;
    r(particleIndices,:) = r(particleIndices,:) + v(particleIndices,:)*tau;
    
    r(particleIndices,1) = mod(r(particleIndices,1)+1000*L,L); %Periodic boundary conditions
    
    [r,v] = collideWithWalls(r,cells,v,R,stddev,rold,p,numParticles,tau);
end

function [r,v] = collideWithWalls(r,cells,v,R,stddev,rold,p,numParticles,tau)
    cj = @(l) ceil((mod(l-1,cells*cells)+1)/cells);
    an = @(x,y) mod((atan2(y,x) + 4*pi),2*pi); % Calculate angle
    transform = @(vector,angle) ([cos(angle) -sin(angle); sin(angle) cos(angle)]*vector')';
    
    
    specularWalls = false;
    
    factor = sqrt(2)*stddev;
    
    for ipart=1:numParticles
       particleIndex = p(ipart); %Real particle index
       if(norm(r(particleIndex,2:3)) > R)
          vz0 = v(particleIndex,3);
          vy0 = v(particleIndex,2);
          % Collided with wall
          angle = an(r(particleIndex,3),r(particleIndex,2));
          
          % Calculate time of collision
          a = vz0^2+vy0^2;
          b = 2*(r(particleIndex,3)*vz0+r(particleIndex,2)*vy0);
          c = r(particleIndex,3)^2 + r(particleIndex,2)^2 - R^2;

          t0 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
          t1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

          dt = max(t0,t1); % Positive root
          
          r(particleIndex,:) = rold(particleIndex,:); % Go back in time
          r(particleIndex,:) = r(particleIndex,:) + v(particleIndex,:)*dt; % Go to wall
          
          if(specularWalls)
              % Inverse the normal velocity component
              vRotated = transform(v(particleIndex,[3,2]),-angle);
              
              % Add sucky friction in parallell velocities
              
              vRotated(1) = -vRotated(1); % z-component
              vRotated(2) = 0.2*vRotated(2); % y-component

              vzy = transform(vRotated,angle);
              vx = 0.2*v(particleIndex,1); % We keep this component as it is

              v(particleIndex,:) = [vx vzy([2,1])];
          else
              % [z,y] are rotated
              vz = -sqrt(-log(rand())) * factor;
              vy = stddev*normrnd(0,1);
              vx = stddev*normrnd(0,1);
              vzy = transform([vz,vy],angle); % Rotate back to normal coordinates
              v(particleIndex,:) = [vx vzy([2,1])];
          end
          
          r(particleIndex,:) = r(particleIndex,:) + v(particleIndex,:)*(tau - dt);
       end
    end
end