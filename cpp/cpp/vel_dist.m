function vel_dist
N = 100000;
r1 = rand(N,1);
r2 = normrnd(0,1,N,1);

v_norm = sqrt(-6/2*3*log(r1));
v_tan = sqrt(3/2*3)*r2;
vsq = v_norm.^2 + v_tan.^2;


end