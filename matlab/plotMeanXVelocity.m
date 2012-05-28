function plotMeanXVelocity(r,v,L,particles)
   res = 1000;
           
   xvel = zeros(res,1);
   count = zeros(res,1);
   for ipart=1:particles
      verCell = ceil(res*r(ipart,2)/L);
      if(verCell > res) verCell = res; end
      xvel(verCell) = xvel(verCell) + v(ipart,1);
      count(verCell) = count(verCell)+1;
   end

   figure(2);
   xvel = smooth(xvel./count);
   
   y = linspace(0,L,res);
   plot(y,xvel)
   title('Mean velocity across the canal')
   xlabel('y')
   ylabel('Mean velocity')
end