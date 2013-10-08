function v =  velocity_profile_plane(file)
    tic
    % hold on
    figure
    data = 0;
    
    if any(strcmp(who,'file'));
        data = dlmread(file);
    else
        data = dlmread('../base_code/statistics/velocity.txt');
    end
    
    
    N_bins_total = size(data,2);
    N_bins_per_dimension = sqrt(N_bins_total);
    N_time_steps = size(data,1);
    time_average = data;
    
    if(N_time_steps > 1)
        time_average = sum(data);
        time_average = time_average/N_time_steps;
    end
    spacing = linspace(0,1,N_bins_per_dimension);
    [x,y] = meshgrid(spacing,spacing);
    % v = zeros(N_bins_per_dimension,N_bins_per_dimension);
    v = reshape(time_average,N_bins_per_dimension,N_bins_per_dimension);
    
%     for k=1:N_bins_total
%         i = floor( (k-1) / N_bins_per_dimension)+1;
%         j = floor(mod(k-1,N_bins_per_dimension)+1);
%         
%         v(i,j) = time_average(k);
%     end
    surf(x,y,v)
    %caxis([0 80])
    toc
end

