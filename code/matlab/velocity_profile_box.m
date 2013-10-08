function velocity_profile_box(file)
    % hold on
    figure
    data = 0;
    
    if any(strcmp(who,'file'));
        data = dlmread(file);
    else
        data = dlmread('../base_code/statistics/velocity.txt');
    end
    
    
    
    N_bins = size(data,2);
    N_time_steps = size(data,1);
    
    velocities = sum(data);
    velocities = velocities/N_time_steps;
    x = linspace(0,1e-6,N_bins);
    j = find(velocities>0);
    x = x(j);
    velocities = velocities(j)
    
    color = [rand() rand() rand()];
    plot(x,velocities,'color',color)
    xlabel('r');
    ylabel('v(r)');
end

