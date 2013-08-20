function plot_velocity_profiles_across_channels(file,data)
    Lz = 1;
    Ly = 0.2;
    
    close all;
    A = 0;
    
    if any(strcmp(who,'data'));
        A = data;
    else
        if any(strcmp(who,'file'));
            A = dlmread(file);
        else
            A = dlmread('../vel_prof/statistics/velocity.txt');
        end
    end
    
    num_timesteps = size(A,1);
    num_bins = size(A,2);
    num_bins_per_dim = sqrt(num_bins);
    length_per_bin_z = Lz/num_bins_per_dim;
    length_per_bin_y = Ly/num_bins_per_dim;
    
    % num_bins_per_dim = num_bins
    
    % Take average over all timesteps
    A = sum(A) / num_timesteps;
    
    z0 = 1+0*num_bins_per_dim;
    num_curves = 10;
    y = zeros(num_curves,num_bins_per_dim);
    
    z_skip = ceil(num_bins_per_dim / num_curves);
    
    for i=1:num_curves
        idx = i-1;
        
        z_index = idx*z_skip*num_bins_per_dim+z0;
        
        y(i,:) = A(z_index:z_index+num_bins_per_dim-1);
        x = linspace(0,Ly,num_bins_per_dim);
        
        figure
        plot(x,y(i,:))
        ylim([0 500])
        xlabel('y [\mu m]');
        ylabel('v(y)');
        title_string = sprintf('Velocity profile at z=%f',z_index/num_bins*Lz);
        title(title_string)
        
    end
    
    % int index = bin_x*num_bins_per_dimension + bin_y;
    
    
end