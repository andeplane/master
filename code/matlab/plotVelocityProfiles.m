function plotVelocityProfiles()
    subplot(1,2,1);
    files = ['/projects/diff1/100/statistics/velocity.txt  ';
             '/projects/diff1/250/statistics/velocity.txt  ';
             '/projects/diff1/500/statistics/velocity.txt  ';
             '/projects/diff1/1000/statistics/velocity.txt ';
             '/projects/diff1/2000/statistics/velocity.txt ';
             '/projects/diff1/5000/statistics/velocity.txt ';
             '/projects/diff1/10000/statistics/velocity.txt'
             ];
    
    for i=1:size(files,1)
        str = strtrim(sprintf('%s',files(i,:)))
        
        velocity_profile_2(str)
    end
    legend('100', '250','500','1000','2000','5000','10000');
    
    subplot(1,2,2);
    files = ['/projects/diff2/100/statistics/velocity.txt  ';
             '/projects/diff2/250/statistics/velocity.txt  ';
             '/projects/diff2/500/statistics/velocity.txt  ';
             '/projects/diff2/1000/statistics/velocity.txt ';
             '/projects/diff2/2000/statistics/velocity.txt ';
             '/projects/diff2/5000/statistics/velocity.txt ';
             '/projects/diff2/10000/statistics/velocity.txt'
             ];
    
    for i=1:size(files,1)
        str = strtrim(sprintf('%s',files(i,:)))
        
        velocity_profile_2(str)
    end
    legend('100', '250','500','1000','2000','5000','10000');
end