function reconstruction = create_approx_solution(W,given_sol,N_proj,d_set,z)

x = given_sol;
n = length(x);
% X = zeros(n,n);
% idx_track = zeros(n,1);
G = W;
Ix = [];
reconstruction = zeros(n,1);
if size(d_set,1)<size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

[reconstruction x Ix] = round_numerically_eq(reconstruction,x,d_set,Ix);
G(:,Ix) = [];

for i_aux = 1:n
    Ix_c = setdiff([1:n],Ix); length(Ix_c)
    rows = sum(G,2) >= N_proj;
    G(~rows,:)=[]; % remove rows of G
    
    if sum(rows) == 0
        reconstruction = round2d_set(reconstruction, given_sol, d_set, Ix_c);
        return
    else
        X0 = rand(size(G,2),1);
        right_h_side = -G*X0;
        [Beta, ~, ~] = cgls_W(G, right_h_side, [], 100, 1e-10);
        Beta = Beta+X0;
        
        if norm(Beta - zeros(length(Beta),1),1)<1e-13
            reconstruction = round2d_set(reconstruction, given_sol, d_set, Ix_c);
            return
        end
        
        Beta_aux = zeros(n,1);
        Beta_aux(Ix_c) = Beta;
        Beta = Beta_aux;
        M_beta = diag(1./Beta);
        t_min = Inf;
        
        for i = 1:length(d_set)
            t_vec = M_beta*(d_set(i)*ones(length(x),1) - x);
            %             t_vec_min = min(abs(t_vec(find(t_vec))));
            t_vec_min = min(t_vec(find(t_vec>0)));
            
            %             if t_vec_min < abs(t_min)
            if t_vec_min < t_min
%                 [index] = find(t_vec==t_vec_min,1);
                %                 if isempty(index)
                %                     [index] = find(t_vec==-t_vec_min,1);
                %                 end
%                 t_min = t_vec(index);
                t_min = t_vec_min;
                t_vec_chosen = t_vec;
            end
        end
        [idx] = find(t_vec_chosen == t_min); % pixels to be removed
        
        x = x + t_min*Beta; % pixel update with rounding error
        x = round2d_set(x, x, d_set, idx); % round to available grey values
    end
    reconstruction(idx) = x(idx); % final pixel value
    Ix = union(Ix,idx); % pixels already defined
    %     [reconstruction x Ix] = round_numerically_eq(reconstruction,x,d_set,Ix);
    G = W;
    G(:,Ix) = []; % pixels removed from the linear system
    %     X(:,i_aux) = x;
    %     idx_track(i_aux,1) = idx;
end

function [reconstruction x Ix] = round_numerically_eq(reconstruction,x,d_set,Ix)
%
for i = setdiff([1:length(x)],Ix)
    for j = 1:length(d_set)
        if abs(x(i)-d_set(j)) <= 1e-14%1/sqrt(length(x))%1e-10
            reconstruction(i) = d_set(j);
            x(i) = d_set(j);
            Ix = union(Ix,i);
        end
    end
end