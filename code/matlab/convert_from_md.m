function convert_from_md(file)
% convert_from_md('~/Dropbox/TightRockModels/Anders/rho05-q3at300-0009.mat');

    M = load(file);
    A = M.rmin>5;
    N = size(A);
    
    fid = fopen('/projects/master/code/worlds/md_m.bin', 'w');
    fwrite(fid, N, 'unsigned char');
    fwrite(fid, A, 'unsigned char');
    fclose(fid);

end

