%MAIN_CREATE_APPROX_RECONSTRUCTION sets the discrete reconstruction problem
%and obtains an approximate solution with bounded projection error.
%Details:
%   Approximate Discrete Reconstruction Algoruthm
%   K.J. Batenburg, W. Fortes, R. Tijdeman
%   Fundamenta Informaticae. Vol 125(3-4), 239-259, 2013.
%
% Wagner Fortes 2014/2015 wfortes@gmail.com

% ------------- parameters:
img_index = 3;  % select image index
img_sz = 64; % select image dimension 32, 64, 128, 256, 512
N_proj_set = [2, 4, 8, 16]; %Set of number of projection angles from 1 to 200
threshold = []; % 1/(img_sz);
loadsol = 1; % boolean: load saved solution in range of d_set
% --------------
% create directions for projection matrix
[dir_a,dir_b]=mkdirvecs(20);

% loads image
P = img_read(img_sz, img_index);
P = reshape(P,img_sz^2,1);
P = double(P);
P = P/norm(P,inf); % only for binary images

if img_index == 7
    d_set = [0;0.50196;1];
else
    d_set = [0;1];
end

for N_proj = N_proj_set
    % Projection matrix
    M = mkmatrix(img_sz,img_sz,dir_a(1:N_proj),dir_b(1:N_proj));
    % Projetion of image P
    Q = M*P;
    
    if loadsol
        filename = strcat('solind_set','Im',img,'-sz',sz,'-proj',proj);
        load(filename);
    else
        x = box_constraint_solver(M, Q, 1e+5*length(Q), 0.1, 'ART', [0;1]);
    end
    %
    %--------------------------------------------------------------------------
    %
    reconstruction = create_approx_solution(M,x,N_proj,d_set,threshold);
    %
    %--------------------------------------------------------------------------
    %
    figure;
    reconstruction = reshape(reconstruction,img_sz,img_sz);
    imshow(reconstruction)
end