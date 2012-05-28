function [r,v] = collideWithEnv(r,sd,cells,v,tau,L)
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells);
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*cells + i;
    friction = 0.1;
    %maxJ = 0;
    for jcell=1:cells*cells
       j = cj(jcell);
       %if j>maxJ
       %    maxJ = j;
       %    sprintf('increased max j: %d',j)
       %end
       
       if j == 1 || j == cells % close to the walls
           particlesInCell = sd(jcell,1);
           
           % Loop through all particles and see if they are colliding
           for ipart = 1:particlesInCell; 
               
              ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index
              
              if r(ip1,2) < 0
                 oldR = r(ip1,:) - v(ip1,:)*tau; %Go back in time
                 dt = oldR(2) / v(ip1,2); % Calculate time until collision
                 
                 r(ip1,:) = oldR + v(ip1,:)*dt; % Move the particles less
                 v(ip1,2) = -v(ip1,2); % Invert velocity
                 v(ip1,:) = v(ip1,:)*friction;
                 % Move the particles the rest of the time
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*(tau - dt);
              end
              
              if r(ip1,2) > L
                 oldR = r(ip1,:) - v(ip1,:)*tau; %Go back in time
                 dt = (L - oldR(2)) / v(ip1,2); % Calculate time until collision
                 
                 r(ip1,:) = oldR + v(ip1,:)*dt; % Move the particles less
                 v(ip1,2) = -v(ip1,2); % Invert velocity
                 v(ip1,:) = v(ip1,:)*friction;
                 % Move the particles the rest of the time
                 r(ip1,:) = r(ip1,:) + v(ip1,:)*(tau - dt);
              end
              
           end
       end
    end
end