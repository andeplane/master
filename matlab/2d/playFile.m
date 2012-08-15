function playFile(filename,timesteps,particles)
    format long;
    L = 1e-6;
    
    if(animate)
        fig1=figure(1);
        winsize = get(fig1,'Position');
        winsize(1:2) = [0 0];
        A=moviein(timesteps,fig1,winsize);
        set(fig1,'NextPlot','replacechildren');
    end
    
    xmax = L*1.05;
    xmin = L-xmax;
    
    fid = fopen(filename,'r');

    tline = fgetl(fid);
    points = zeros(1,3);
    numPoints = 0;
    time = 0;
    r = zeros(timesteps,3);
    
    while ischar(tline) 
        if(numel(tline) > 0 && strcmp(tline(1),'a'))
            
            % r(numPoints,:) = str2num(tline(2:length(tline)));
        end

        tline = fgetl(fid);
    end
    fclose(fid);
    
    
    for i=1:timesteps
       if(energyVsTime) E(i) = 0.5*sum(sum(v.^2,2)); end
        
       if(animate)
           % Clear figure and plot this timestep
           figure(1);
           clf(fig1);
           plot(r(:,1),r(:,2),'d');
           hold on;
           plot([0 L],[L,L],'r');
           plot([0 L],[0,0],'r');
           plot([L L],[0,L],'c');
           plot([0 0],[0,L],'c');
           xlabel('x [m]');
           ylabel('y [m]');
           axis('equal')
           % axis([xmin xmax xmin xmax]);

           title(sprintf('t = %f ns',10^9*i*tau));
           A(:,i)=getframe(fig1,winsize); % save to movie
       end
       v_before = v;
       
       fprintf(fid,'%d\n a\n',particles);
       fprintf(fid,'H %f %f 0\n',r(:,1:2)*10^6);
       
       sortData = sortParticles(r,L,cells,particles,sortData); %Put particles in cells
       [r,v] = mover(r,sortData,cells,v,tau,L,v0); % Move and collide with walls
       [col, v, vrmax, selxtra] = collide(v,vrmax,selxtra,coeff,sortData,cells,L); % Collide with other particles
       deltaV = v - v_before;
       
       % [P,Pmax] = calculatePressure(sortData,cells,v,mass,L,plotPressure,Pmax,eff_num,deltaV);
       
       
       coltot = coltot + col; % Increase total collisions
       
       if mod(i,100) < 1
           dCol = coltot - oldColTot;
           oldColTot = coltot;
           if meanXVelocity plotMeanXVelocity(r,v,L,particles); end
           
           
           sprintf('Done %d of %d steps. %d collisions (%d new)',i,timesteps,coltot,dCol)
       end
    end
    fclose(fid);
    sprintf('Done!')
    
    if(energyVsTime)
        t = linspace(0,timesteps*tau,timesteps);
        plotEnergy(t,E);
    end
    
end