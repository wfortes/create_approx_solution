function reconstruction = create_approx_solution(W, given_sol, N_proj, d_set, thresh)
%CREATE_APPROX_SOLUTION obtains an approximate solution with bounded 
%projection error.
%   Details:
%   Approximate Discrete Reconstruction Algoruthm
%   K.J. Batenburg, W. Fortes, R. Tijdeman
%   Fundamenta Informaticae. Vol 125(3-4), 239-259, 2013.
%
%   Input:
%   W matrix of the reconstruction problem
%   GIVEN_SOL requires a initial solution on the range of the D_SET
%   N_PROJ number of projection angles
%   D_SET set of discrete values of the original image
%   THRESH threshold for rounding approx solution at each iteration
%
%   Output:
%   RECONSTRUCTION is the reconstructed discrete image with projection
%   guaratees
%     
%
% Wagner Fortes 2014/2015 wfortes@gmail.com

x = given_sol;
n = length(x);
A = W;
Ix = [];
iterations = 0;
reconstruction = zeros(n,1);

if size(d_set,1) < size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

[reconstruction, x, Ix] = round_numerically_eq(reconstruction,x,d_set,Ix,thresh);
A(:,Ix) = [];

Ix_c = setdiff([1:n],Ix); 
% length(Ix_c)
rows = sum(A,2) < N_proj;
A(rows,:)=[]; % remove rows of A
rows0 = sum(A,2) == 0; % identify the null rows
A(rows0,:)=[]; % remove rows of A

while ~isempty(A)
    
    iterations = iterations + 1;
    X0 = rand(size(A,2),1);
    right_h_side = -A*X0;
    Beta = ls_solver(A, right_h_side, 'cgls_W', [], []);
    Beta = Beta+X0;
    
    if norm(Beta - zeros(length(Beta),1),1) < 1e-13
        break
    end
    
    Beta_aux = zeros(n,1);
    Beta_aux(Ix_c) = Beta;
    Beta = Beta_aux;
    M_beta = diag(1./Beta);
    t_min = Inf;
    t_vec_chosen =[];
    for i = 1:length(d_set)
        t_vec = M_beta*(d_set(i)*ones(length(x),1) - x);
        t_vec_min = min(t_vec(t_vec>0));
        if t_vec_min < t_min
            t_min = t_vec_min;
            t_vec_chosen = t_vec;
        end
    end
    if isempty(t_vec_chosen)
        error('shouldnt be empty')
    end
    [idx] = find(t_vec_chosen == t_min); % pixels to be removed
    
    x = x + t_min*Beta; % pixel update with rounding error
    x = round2d_set(x, x, d_set, idx); % round to available grey values
    
%     reconstruction(idx) = x(idx); % final pixel value
    if isempty(thresh)
        reconstruction = x; % final pixel value
    else
        [reconstruction, x, Ix] = round_numerically_eq(reconstruction,x,d_set,Ix,thresh);
    end
    
    Ix = union(Ix,idx); % pixels already defined
    
    A = W;
    A(:,Ix) = []; % pixels removed from the linear system
    Ix_c = setdiff([1:n],Ix); 
    
    if mod(iterations,50) == 0
    fprintf('%d pixels left\n', length(Ix_c));
    end
    
    rows = sum(A,2) < N_proj;
    
    A(rows,:)=[]; % remove rows of A
    rows0 = sum(A,2) == 0; % identify the null rows
    A(rows0,:)=[]; % remove rows of A
end
reconstruction = round2d_set(reconstruction, given_sol, d_set, Ix_c);

function [reconstruction, x, Ix] = round_numerically_eq(reconstruction,x,d_set,Ix,thresh)
%
if isempty(thresh)
    thresh = 1e-14;
end
for i = setdiff([1:length(x)],Ix)
    for j = 1:length(d_set)
        if abs(x(i)-d_set(j)) <= thresh
            reconstruction(i) = d_set(j);
            x(i) = d_set(j);
            Ix = union(Ix,i);
        end
    end
end