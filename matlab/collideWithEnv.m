function [r,v] = mover(r,sd,cells,v,tau,L,mpv)
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells);
    cj = @(k) ceil(k/cells);
    %cell number as function of i, j
    ck = @(i,j) (j-1)*cells + i;
    friction = 0.1;
    
    direction = [1, -1];
    ywall = [0, L];
    stddev = mpv/sqrt(2);
    
    yold = r;
              
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
              flag = 0;
              
              if(r(ip1,2) <= 0) flag = 1; end
              if(r(ip1,2) >= L) flag = 2; end
              
              if(flag > 0)
                % we collided
                v(ip1,2) = direction(flag)*sqrt(-log(1-rand())) * mpv;
                v(ip1,1) = stddev*normrnd(0,1);
                sprintf('New y=%f, old y=%f, diff=%f',
                dtr = tau*(r(ip1,2)-ywall(flag))/(r(ip1,2)-yold(ip1,2));
                
                r(ip1,2) = ywall(flag) + v(ip1,2)*dtr;
                sprintf('We collided (flag=%f). dtr=%f',flag,dtr)
              end
           end
       end
    end
end