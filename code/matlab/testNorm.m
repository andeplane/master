function testNorm()
    tx = @(v) v(2);
    ty = @(v) v(1);
    tanvec = @(v) [tx(v) -ty(v)];

    %r = rand(3,3);
    %m = r>0.5
    m = [ 1 0 0;
          1 1 0;
          0 0 0
        ]
     
    v = [0 0];
    
    for i=1:3
        for j=1:3
           k = i-2;
           l = j-2;
           
           v(1) = v(1)-m(i,j)*k;
           v(2) = v(2)+m(i,j)*l;
        end
    end
    v = v/norm(v)
    t = tanvec(v)
    
    figure
    subplot(2,2,1)
    imshow(m)
    subplot(2,2,2)
    x0 = [0 0];
    vectarrow(x0,v);
    hold on
    %vectarrow(x0,t);
    axis([-1 1 -1 1])
    %arrow([0,v(1)],[0,v(2)]);
end