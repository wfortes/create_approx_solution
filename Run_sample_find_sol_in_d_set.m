%RUN_SAMPLE_FIND_SOL_IN_D_SET sets the discrete tomography problem and
%finds a solution that belongs to a range given by the extreme values of
%given set. Save solution on file.
%
% Wagner Fortes 2014/2015 wfortes@gmail.com

img_index = 7;
img_sz = 64;
N_proj_set = [4, 8, 16];

[dir_a,dir_b]=mkdirvecs(20);

P = img_read(img_sz, img_index);
P = reshape(P,img_sz^2,1);
P = double(P);
P = P/norm(P,inf); % only for binary images


for N_proj = N_proj_set
    
    M = mkmatrix(img_sz,img_sz,dir_a(1:N_proj),dir_b(1:N_proj));
    Q = M*P;
    
    x = box_constraint_solver(M, Q, 1e+5*length(Q), 0.1, 'ART', [0;1]);
    
    img = num2str(img_index);
    sz = num2str(img_sz);
    proj = num2str(N_proj);
    
    filename = strcat('solind_set','Im',img,'-sz',sz,'-proj',proj);
    save(filename,'x');
end

