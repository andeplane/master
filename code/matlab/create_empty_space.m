N = [10 10 10];
A = zeros(N);

fid = fopen('/projects/master/code/worlds/empty_m.bin', 'w');
fwrite(fid, N, 'unsigned char');
fwrite(fid, A, 'unsigned char');
fclose(fid);