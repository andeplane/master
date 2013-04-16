function plotVelocityProfiles()
    base = '/projects/states/0.4';
    
    subplot(1,2,1);
    density = 1e25;
    atoms_per_molecule = [100, 250, 500, 1000, 2000, 5000, 10000];
    
    for i=1:size(atoms_per_molecule,2)
        file = sprintf('%s/%e-%d/statistics/velocity.txt',base,density,atoms_per_molecule(i) )
        
        velocity_profile_2(file)
    end
    legend('100', '500','1000','2000','5000','10000');
    xlabel('y')
    ylabel('v(y)')
    title('rho=1e25/m^3')
    
    subplot(1,2,2);
    density = 1e26;
    atoms_per_molecule = [100, 500, 1000, 2000, 5000, 10000];
    
    for i=1:size(atoms_per_molecule,2)
        file = sprintf('%s/%e-%d/statistics/velocity.txt',base,density,atoms_per_molecule(i) )
        
        velocity_profile_2(file)
    end
    legend('100', '500','1000','2000','5000','10000');
    xlabel('y')
    ylabel('v(y)')
    title('rho=1e26/m^3')
end