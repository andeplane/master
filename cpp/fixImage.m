function fixImage( file )
%FIXIMAGE Summary of this function goes here
%   Detailed explanation goes here

img = imread(file);
imshow(img)
figure
points = [
123,185
123,185
];

% points = [points(:,2) points(:,1)];
points(:,1) = size(img,1)-points(:,1);

for i=1:length(points)
    img(points(i,1),points(i,2),:) = 0;
    img(points(i,1),points(i,2),1) = 255;
end

idx_i = @(i) mod(i-1 + size(img,2),size(img,2))+1;
idx_j = @(j) mod(j-1 + size(img,1),size(img,1))+1;

bool_img = zeros(size(img,1),size(img,2));

for i=1:size(img,2)
    for j=1:size(img,1)
       bool_img(j,i) = sum(img(j,i,:) > 0);
    end
end

for i=1:size(img,2)
    for j=1:size(img,1)
        hasNormals = false;
        normal = zeros(2,1);
        for k=-1:1
            for l=-1:1
               %i_i = idx_i(i+k);
               %i_j = idx_j(j+l);
               i_i = i+k;
               i_j = j+l;
               if(i_j < 1)           i_j = size(img,1); end
               if(i_j > size(img,1)) i_j = 1; end
               if(i_i < 1)           i_i = size(img,2); end
               if(i_i > size(img,2)) i_i = 1; end
               
               normal(1) = normal(1) - bool_img(j,i)*k;
               normal(2) = normal(2) - bool_img(j,i)*l;
               if(normal'*normal > 0) hasNormals = true; end
            end
        end
        if(bool_img(j,i) && ~hasNormals)
           % img(j,i,:) = 0;
        end
    end
    
    sprintf('%.2f',100*i/size(img,2))
end

imwrite(img,'test2.bmp');
imshow(img)
end