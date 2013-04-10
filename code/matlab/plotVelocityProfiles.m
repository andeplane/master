function plotVelocityProfiles(data)
    %data = dlmread('../compiled/release/velocity.txt');
    size(data)
    num_samples = size(data,1);
    num_velocities_total = size(data,2);
    num_velocities_vectors_total = num_velocities_total/3; % 3 each 
    
    num_cells_per_dim = sqrt(num_velocities_vectors_total);
    v = zeros(num_cells_per_dim,num_cells_per_dim);
    vr = zeros(100,1);
    vr_count = zeros(100,1);
    
    for n=1:num_samples
        for i=1:num_cells_per_dim
            for j=1:num_cells_per_dim
                index = (i-1)*num_cells_per_dim + j-1;
                
                this_v = data(n,3*index+3);
                v(i,j) = v(i,j) + this_v;
                dr = sqrt( (i-26)^2 + (j-26)^2 );
                
                vr_index = round(dr)+1;
                vr(vr_index) = vr(vr_index) + this_v;
                vr_count(vr_index) = vr_count(vr_index) + 1;
            end
        end
    end
    
    for i=1:100
        if(vr_count(i) > 0)
            vr(i) = vr(i) / vr_count(i); 
        end
    end
    x = linspace(0,0.5,length(vr));
    plot(x,vr,'r')
    xlabel('r [\mu m]');
    ylabel('v [m^s]');
    
    figure
    for i=1:num_cells_per_dim
        for j=1:num_cells_per_dim
            v(i,j) = v(i,j)/num_samples;
        end
    end
    
    surf(v);
    xlabel('y');
    ylabel('z');
   
end