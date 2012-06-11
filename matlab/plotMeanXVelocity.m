   function plotMeanXVelocity(r,v,R,numParticles,p)
   
   res = 50;
   
   rvel = zeros(res,1);
   count = zeros(res,1);
   numRandomSamples = numParticles;
   sampleIndices = randsample(p(1:numParticles),numRandomSamples);
   
   for ipart=1:numRandomSamples
      particleIndex = sampleIndices(ipart);
      
      verCell = ceil(res*norm(r(particleIndex,2:3))/R);

      rvel(verCell) = rvel(verCell) + norm(v(particleIndex,:));
      count(verCell) = count(verCell)+1;
   end

   figure(2);
   xvel = smooth(rvel./count);
   
   y = linspace(0,1,res);
   plot(y,xvel)
   title('Mean velocity across the tube')
   xlabel('r/R')
   ylabel('Mean velocity')
end