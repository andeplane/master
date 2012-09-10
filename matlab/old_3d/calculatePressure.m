function [P,Pmax] = calculatePressure(sd,cells,v,mass,L,doPlotPressure,Pmax,eff_num,deltaV)
    %cell number i and j as a function of k
    ci = @(k) mod(k-1,cells)+1;
    cj = @(k) ceil(k/cells);
    
    V = L^2/cells^2;
    
    P = zeros(cells,cells);
    P1 = zeros(cells,cells);
    
    for jcell = 1:cells*cells
        number = sd(jcell,1); %Number of particles in this cell
        
        if number < 1
           continue; 
        end
        
        i = ci(jcell);
        j = cj(jcell);
        
        sumv = 0;
        
        for ipart=1:number
           ip1 = sd(ipart+sd(jcell,2)-1,3); % Actual particle index
           sumv = sumv+norm(v(ip1,:))^2;
           % P1(i,j) = P1(i,j) + norm(deltaV(ip1,:));
        end
        
        P1(i,j) = number;
        
        P(i,j) = sumv*eff_num*mass/(L*L);
    end
    
    % P = P1*mass/(L^3);
    
    newMax = max(P(:));
    if(newMax > Pmax)
       Pmax = newMax; 
    end
    
    if(doPlotPressure)
        plotPressure(P,L,Pmax);
    end
end

function plotPressure(P,L,Pmax)
    
    figure(4);
    x = 0:L;
    y = 0:L;
    imagesc(x,y,P');
    colorbar;
    caxis([0 Pmax])
    sumP = sum(P(:))
end