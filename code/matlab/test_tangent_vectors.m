function test_tangent_vectors()
    N = 10000;
    
    normal = [0 1 0];
    tangent_vectors_1 = rand(N,3);
    tangent_vectors_2 = zeros(N,3);
    gaussian_numbers = normrnd(0,1,N,2);
    skewed_gaussian_numbers = sqrt(-2*log(rand(N,1)));
    velocities = zeros(N,3);
    v_norms = zeros(N,1);
    
    for i=1:N
        tangent_vectors_1(i,:) = tangent_vectors_1(i,:) / norm(tangent_vectors_1(i,:),2);
        tangent_vectors_2(i,:) = cross(normal,tangent_vectors_1(i,:));
        
        velocities(i,1) = normal(1)*skewed_gaussian_numbers(i) + tangent_vectors_1(i,1)*gaussian_numbers(i,1) + tangent_vectors_2(i,1)*gaussian_numbers(i,2);
        velocities(i,2) = normal(2)*skewed_gaussian_numbers(i) + tangent_vectors_1(i,2)*gaussian_numbers(i,1) + tangent_vectors_2(i,2)*gaussian_numbers(i,2);
        velocities(i,3) = normal(1)*skewed_gaussian_numbers(i) + tangent_vectors_1(i,3)*gaussian_numbers(i,1) + tangent_vectors_2(i,3)*gaussian_numbers(i,2);
        v_norms(i) = norm(velocities(i,:),2);
    end
    
    hist(v_norms,30)
    
    
    
end