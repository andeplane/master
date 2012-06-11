   function plotMeanXVelocity(r,v,R,numParticles,p)
   
   res = 1000;
   
   rvel = zeros(res,1);
   count = zeros(res,1);
   numRandomSamples = 2000;
   sampleIndices = randsample(p(1:numParticles),numRandomSamples);
   
   for ipart=1:numRandomSamples
      particleIndex = sampleIndices(ipart);
      
      if(abs(r(particleIndex,2) < R/10))
          %verCell = ceil(res*norm(r(particleIndex,2:3))/R);
          verCell = ceil(res*(R+r(particleIndex,3))/(2*R));
          
          rvel(verCell) = rvel(verCell) + v(particleIndex,1);
          count(verCell) = count(verCell)+1;
      end
      
      
   end

   figure(2);
   xvel = smooth(rvel./count);
   
   y = linspace(0,1,res);
   plot(y,xvel)
   title('Mean velocity across the tube')
   xlabel('r/R')
   ylabel('Mean velocity')
end