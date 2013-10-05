N = [100 100 100];
A = zeros(N);

fid = fopen('/projects/master/code/worlds/empty_m.bin', 'w');
fwrite(fid, N, 'uint');
fwrite(fid, A, 'unsigned char');
fclose(fid);