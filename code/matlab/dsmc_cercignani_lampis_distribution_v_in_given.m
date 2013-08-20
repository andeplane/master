function [v_out, f] = dsmc_cercignani_lampis_distribution_v_in_given(v_in, alpha_n)
    v_out_max = 1500;
    k = 1.3806488e-23;
    T = 300;
    m = 6.6e-26;
    factor = sqrt(2*k*T/m);
    v_in = v_in/factor;
    v_out = linspace(0,v_out_max/factor,100);
    
    bessel_argument = 2/alpha_n*sqrt(1-alpha_n)*v_in.*v_out;
    
    C = 2/alpha_n;
    a = exp(-(v_out - sqrt(1-alpha_n)*v_in).^2/alpha_n);
    b = v_out.*exp(-bessel_argument).*besseli(0,bessel_argument);
    f = C*a.*b;
    %Z = 2*v_out/alpha_n.*exp( -(v_out.^2 + (1-alpha_n)*v_in.^2)/alpha_n).*besseli(0,bessel_argument);
    
    v_out = v_out*factor;
    plot(v_out, f);
    hold all;
    plot(v_out, a);
    plot(v_out, b);
    legend('f(x)','a(x)','b(x)')
    
    xlabel('v_{out}')
    ylabel('P')
end