mfp = [0.0245,0.0573,0.0700,0.0859,0.1718, 0.3435];
close all;
for i=1:length(mfp)
   velocity_profile_2(mfp(i));
end

legend('Kn=0.03','Kn=0.07','Kn=0.9','Kn=0.11','Kn=0.22','Kn=0.43')