function reconstruction = create_approx_solution_piecewise(W,given_sol,N_proj,d_set,piece)

% this is only implemented for square images and divisible by the piece

n = length(given_sol);
mm = n/piece^2;
m = sqrt(mm);
reconstruction_partial = given_sol;

if size(d_set,1)<size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

for i_pieces_columns = 1:m
    for i_pieces_rows = 1:m
        Ix = [];
        pixels_in_piece = [];
        
        for i_row = 0:piece-1
            pixels_in_piece_aux = i_row*piece*m+1+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1):i_row*piece*m+piece+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1);
            pixels_in_piece = union(pixels_in_piece,pixels_in_piece_aux);
        end
        x = given_sol(pixels_in_piece);
        G = W(:,pixels_in_piece);
        
        [x Ix] = round_numerically_eq(x,d_set,Ix);
        G(:,Ix) = [];
        Ix_c = setdiff(1:piece^2,Ix); length(Ix_c)
        rows = sum(G,2) >= N_proj;
        
        while sum(full(rows)) ~= 0
            
            G(~rows,:)=[]; % remove rows of G
            
            X0 = rand(size(G,2),1);
            right_h_side = -G*X0;
            [Beta, ~, ~] = cgls_W(G, right_h_side, [], 100, 1e-15);
            Beta = Beta+X0;
            if norm(Beta - zeros(length(Beta),1),1)<1e-14
                reconstruction_partial(pixels_in_piece) = x;
                break
            end
%             if norm(Beta-zeros(length(Beta),1),1)<1e-13
%                 error('Beta errado')
%             end
            
            Beta_aux = zeros(piece^2,1);
            Beta_aux(Ix_c) = Beta;
            Beta = Beta_aux;
            M_beta = pinv(diag(Beta));
            t_min = Inf;
            
            for i = 1:length(d_set)
                t_vec = M_beta*(d_set(i)*ones(length(x),1) - x);
                t_vec_min = min(t_vec(find(t_vec>0)));
                
                if t_vec_min < t_min
                    t_min = t_vec_min;
                    t_vec_chosen = t_vec;
                end
            end
            [idx] = find(t_vec_chosen == t_min); % pixels to be removed
            x = x + t_min*Beta; % pixel update with rounding error
            x = round2d_set(x, x, d_set, idx); % round to available grey values
            
%             reconstruction_partial(pixels_in_piece) = x;
            Ix = union(Ix,idx); % pixels already defined
            G = W(:,pixels_in_piece);
            G(:,Ix) = []; % pixels removed from the linear system
            Ix_c = setdiff(1:piece^2,Ix); length(Ix_c)
            rows = sum(G,2) >= N_proj;
            
%             reconstruction_partial(pixels_in_piece) = x;
%             norm(W(:,pixels_in_piece)*(x-given_sol(pixels_in_piece)),inf)
%             norm(W*(reconstruction_partial-given_sol),inf)
            if norm(W(:,pixels_in_piece)*(x-given_sol(pixels_in_piece)),inf)>(N_proj-1)*(d_set(end)-d_set(1))+1e-15
                error('Esta errado')
            end
        end
        reconstruction_partial(pixels_in_piece) = x;
    end
end

reconstruction = create_approx_solution_aux(W,given_sol,reconstruction_partial,N_proj,d_set);

function reconstruction = create_approx_solution_aux(W,given_sol,intermediate,N_proj,d_set)

x = intermediate;
n = length(x);
G = W;
Ix = [];
if size(d_set,1)<size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

[x Ix] = round_numerically_eq(x,d_set,Ix);
G(:,Ix) = [];
Ix_c = setdiff([1:n],Ix); length(Ix_c)
rows = sum(G,2) >= N_proj;

while sum(rows) ~= 0
    G(~rows,:)=[]; % remove rows of G
    
    X0 = rand(size(G,2),1);
    right_h_side = -G*X0;
    [Beta, ~, ~] = cgls_W(G, right_h_side, [], 100, 1e-15);
    Beta = Beta+X0;
    
    if norm(Beta - zeros(length(Beta),1),1)<1e-14
        reconstruction = round2d_set(x, given_sol, d_set, Ix_c);
        return
    end
    
    Beta_aux = zeros(n,1);
    Beta_aux(Ix_c) = Beta;
    Beta = Beta_aux;
    M_beta = pinv(diag(Beta));
    t_min = Inf;
    
    for i = 1:length(d_set)
        t_vec = M_beta*(d_set(i)*ones(length(x),1) - x);
        t_vec_min = min(t_vec(find(t_vec>0)));
        
        if t_vec_min < t_min
            t_min = t_vec_min;
            t_vec_chosen = t_vec;
        end
    end
    [idx] = find(t_vec_chosen == t_min); % pixels to be removed
    
    x = x + t_min*Beta; % pixel update with rounding error
    x = round2d_set(x, x, d_set, idx); % round to available grey values
    
    Ix = union(Ix,idx); % pixels already defined
    G = W;
    G(:,Ix) = []; % pixels removed from the linear system
    Ix_c = setdiff([1:n],Ix); length(Ix_c)
    rows = sum(G,2) >= N_proj;
    
end
reconstruction = round2d_set(x, given_sol, d_set, Ix_c);

function [x Ix] = round_numerically_eq(x,d_set,Ix)
%
for i = setdiff([1:length(x)],Ix)
    for j = 1:length(d_set)
        if abs(x(i)-d_set(j)) <= 1e-14%1/sqrt(length(x))
            x(i) = d_set(j);
            Ix = union(Ix,i);
        end
    end
end