function testNorm()
    %r = rand(3,3);
    %m = r>0.5
    m = [     0        0   1.0000;
            1.0000   1.0000   1.0000;
            1.0000   1.0000   1.0000]
     
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
    
    figure
    subplot(2,2,1)
    imshow(m)
    subplot(2,2,2)
    x0 = [0 0];
    vectarrow(x0,v);
    axis([-1 1 -1 1])
    %arrow([0,v(1)],[0,v(2)]);
end