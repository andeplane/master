function density()
    figure

    density = dlmread('../base_code/statistics/density.txt');
    
    num_bins_total = length(density);
    num_bins_per_dimension = sqrt(num_bins_total);
    
    density = reshape(density,num_bins_per_dimension,num_bins_per_dimension);
    
    x = linspace(0,10,num_bins_per_dimension);
    y = linspace(0,10,num_bins_per_dimension);
    
    [x,y] = meshgrid(x,y);
    size(x)
    size(y)
    size(density)
    surf(x,y,density);
    %caxis([1.5 2.5])
end    