% A sampling technique to sample from the complicated distribution in the
% Cercignani-Lampis model for gas-surface interactions.
% Sampling technique source: Computer Simlation of Liquids - M. Allen,
% Oxford 1989 (p 350)

function dsmc_sample_cercignani_lampis_distribution(N, v_in, alpha_n)
    close all;
    k = 1.3806488e-23;
    T = 300;
    m = 6.6e-26;
    factor = sqrt(2*k*T/m);
    v_in = v_in/factor;
    
    % Step 1: generate a uniform random variate ksi on (0,1)
    uniform_rnd = rand(N,1);
    ksi = uniform_rnd;
    % Step 2: generate zeta randomly on distribution a(x), in this case
    % normal distribution.
    std_dev = sqrt(alpha_n/2);
    mu = sqrt(1-alpha_n)*v_in;
    gauss_rnd = normrnd(mu,std_dev,N,1);
    v_out = gauss_rnd;
    
    % Step 3: if ksi < b(zeta), zeta is random on p(x)
    bessel_argument = 2/alpha_n*sqrt(1-alpha_n)*v_in.*v_out;
    
    b = v_out.*exp(-bessel_argument).*besseli(0,bessel_argument);
    accepted_indices = b>=ksi;
    v_out = v_out(accepted_indices);
    
    sprintf('%d of %d (%f%ff) success',length(v_out), N, length(v_out)/N)
    
    subplot(3,2,1)
    [num,X] = hist(factor*v_out,30);
    hist(factor*v_out,30);
    subplot(3,2,2)
    [v_out, f] = dsmc_cercignani_lampis_distribution_v_in_given(v_in*factor, alpha_n);
    f = f * max(num)/max(f);
    subplot(3,2,1);
    hold on
    plot(v_out, f,'ro');
    subplot(3,2,3);
    A = dlmread('/projects/master/code/cpp/result.txt');
    [num,X] = hist(A,30);
    hist(A,30);
    hold on
    f = f * max(num)/max(f);
    plot(v_out, f,'ro');
end


function plot_b(alpha_n)
    k = 1.3806488e-23;
    T = 300;
    m = 6.6e-26;
    factor = sqrt(2*k*T/m);
    
    v_in = linspace(-5,0,100);
    v_out = linspace(0,5,100);
    
    [V_IN,V_OUT] = ndgrid(v_in,v_out);
    
    bessel_argument = 2/alpha_n*sqrt(1-alpha_n)*V_IN.*V_OUT;
    
    b = V_OUT.*exp(-bessel_argument).*besseli(0,bessel_argument);
    
    h_surf = surf(factor*V_IN, factor*V_OUT, b);
    xlabel('v_{in}')
    ylabel('v_{out}')
    zlabel('b')
    set(h_surf,'LineStyle','none')
    view(2)
end