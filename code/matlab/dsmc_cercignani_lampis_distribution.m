function dsmc_cercignani_lampis_distribution()
    close all;
    k = 1.3806488e-23;
    T = 100;
    m = 6.6e-26;
    factor = sqrt(2*k*T/m);
    alpha_n = 0.5;
    v_in = linspace(-10,0,100);
    v_out = linspace(0,10,100);
    
    [V_IN,V_OUT] = ndgrid(v_in,v_out);
    
    bessel_argument = 2/alpha_n*sqrt(1-alpha_n)*V_IN.*V_OUT;
    
    Z = 2*V_OUT/alpha_n.*exp( -(V_OUT.^2 + (1-alpha_n)*V_IN.^2)/alpha_n).*besseli(0,bessel_argument);
    
    h_surf = surf(factor*V_IN, factor*V_OUT, Z);
    xlabel('v_{in} [m/s]', 'fontsize', 15);
    ylabel('v_{out} [m/s]', 'fontsize', 15);
    zlabel('p(v_{in},v_{out})', 'fontsize', 15);
    title('Cercignani-Lampis normal component distribution ', 'fontsize',13);
    set(h_surf,'LineStyle','none');
    view(2)
    colorbar()
end