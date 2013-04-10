function velocity_profile_2(mfp)
    hold on
    A = dlmread('../compiled/release/velocity.txt');
    %kn = mfp/0.8
    %file = sprintf('../compiled/release/data/velocity_%.4f.txt',mfp);
    %file
    %A = dlmread(file);
    
    N_bins = size(A,2);
    N_time_steps = size(A,1);
    
    B = sum(A);
    B = B/N_time_steps;
    B = B/max(B);
    x = linspace(0,1,N_bins);
    %j = find(B>0);
    
    color = [rand() rand() rand()];
    plot(x,B,'color',color)
    xlabel('r');
    ylabel('v(r)');
    
    
end

