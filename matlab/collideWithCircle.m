function [r,v] = collideWithCircle(r,sd,cells,v,tau,radius,xc,yc,L,particles)
    for jcell=1:cells*cells
       particlesInCell = sd(jcell,1);

       % Loop through all particles and see if they are colliding
       
       for ipart = 1:particlesInCell; 
       %for ipart = 1:particles 
          ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index
          %ip1 = ipart;

          rFromCircle = r(ip1,:) - [xc,yc,0];
          
          if(norm(r(ip1,:) - [xc,yc,0]) < radius)
             %collided

             oldR = r(ip1,:) - v(ip1,:)*tau; %Go back in time
             % Calculate time of colision
             x0 = oldR(1);
             y0 = oldR(2);
             vx = v(ip1,1);
             vy = v(ip1,2);

             a = v(ip1,1)*v(ip1,1) + v(ip1,2)*v(ip1,2);
             b = 2*(x0*vx + y0*vy - xc*vx - yc*vy);
             c = -radius^2 - 2*xc*x0 - 2*yc*y0 + x0^2 + y0^2 + yc^2 + xc^2;
             t0 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
             t1 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

             t = min(t0,t1);

             r(ip1,:) = oldR + v(ip1,:)*t;
             v(ip1,:) = -v(ip1,:);

             r(ip1,:) = r(ip1,:) + v(ip1,:)*(tau-t);
             
             % It is possible that we send our particle outside the box,
             % So for now, just cap the y-position inside [0,L]
             r(ip1,2) = min(r(ip1,2),L);
             r(ip1,2) = max(r(ip1,2),0);
             
          end
          
       end
    end
    
    for ip1 = 1:particles
       if(norm(r(ip1,:) - [xc,yc,0]) < radius) 
             sprintf('Particle %d is still inside the circle. Distance after: %f',ip1,norm(r(ip1,:) - [xc,yc,0])/radius)
         end 
    end

end