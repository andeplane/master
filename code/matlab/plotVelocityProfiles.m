function plotVelocityProfiles()
    files = ['../compiled/release/data/velocity_200k100k.txt';
             '../compiled/release/data/velocity_150k50k.txt ';];
    
    for i=1:size(files,1)
        str = strtrim(sprintf('%s',files(i,:)))
        
        velocity_profile_2(str)
    end
end