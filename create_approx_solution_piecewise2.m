function reconstruction = create_approx_solution_piecewise2(W,given_sol,N_proj,d_set,piece)

% this is only implemented for square images and divisible by the piece

n = length(given_sol);
mm = n/piece^2;
m = sqrt(mm);
reconstruction_partial = given_sol;
Ix_global = [];
res = [];

if size(d_set,1)<size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

[x Ix] = round_numerically_eq(given_sol,d_set,[]);% global
Ix_global = union(Ix_global,Ix);

for i_pieces_columns = 1:m
    for i_pieces_rows = 1:m
        Ix = [];
        pixels_in_piece = [];
        
        for i_row = 0:piece-1
            pixels_in_piece_aux = i_row*piece*m+1+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1):i_row*piece*m+piece+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1);
            pixels_in_piece = union(pixels_in_piece,pixels_in_piece_aux); %global coordinates
        end
        x = given_sol(pixels_in_piece);
        
%         [x Ix] = round_numerically_eq(x,d_set,Ix); % local
        Ix_c = setdiff(1:piece^2,Ix);
        
        %         Ix_global = union(Ix_global,pixels_in_piece(Ix));
        G = W;
        G(:,Ix_global) = []; % remove columns from W
        rows = sum(G,2) < N_proj;
        
        A = W(:,pixels_in_piece); % remove columns from W
        A(:,Ix) = []; % remove columns local
        
        A(rows,:)=[]; % remove rows of A according to the set G
        rows0 = sum(A,2) == 0; % identify the null rows
        A(rows0,:)=[]; % remove null rows of A
        
        while ~isempty(G)||~isempty(A)
            
            X0 = rand(size(A,2),1);
            right_h_side = -A*X0;
            [Beta, ~, ~] = cgls_W(A, right_h_side, [], 1000, 1e-15);
            Beta = Beta+X0;
            res = res + norm(A*Beta,inf);
            if norm(Beta - zeros(length(Beta),1),1)<1e-14
                warning('null solution')
                break
            end
            
            Beta_aux = zeros(piece^2,1);
            Beta_aux(Ix_c) = Beta; % local
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
            
            reconstruction_partial(pixels_in_piece) = x;
            H=W;
            H(rows,:)=[];
            for au = 1:size(H,1)
                if abs(H(au,:)*(reconstruction_partial-given_sol))>res % MAKE RES MORE PRECISE
                    error('1');
                end
            end
            HH=W;
            HH(~rows,:)=[];
            for au = 1:size(HH,1)
                if abs(HH(au,:)*(reconstruction_partial-given_sol))> (N_proj-1)*(d_set(end)-d_set(1))+res % MAKE RES MORE PRECISE
                    error('2');
                end
            end
            
            Ix = union(Ix,idx); % pixels already defined
            Ix_c = setdiff(1:piece^2,Ix); length(Ix_c)
            
            Ix_global = union(Ix_global,pixels_in_piece(idx));
            G = W;
            G(:,Ix_global) = []; % pixels removed from the set G
            rows = sum(G,2) < N_proj;
            
            A = W(:,pixels_in_piece); % remove columns of A
            A(:,Ix) = []; % remove columns of A (local)
            
            A(rows,:)=[]; % remove rows of A according to the set G
            rows0 = sum(A,2) == 0; % identify the null rows
            A(rows0,:)=[]; % remove null rows of A
        end
        reconstruction_partial(pixels_in_piece) = x;
    end
end

reconstruction = create_approx_solution_aux(W,given_sol,reconstruction_partial,N_proj,d_set,res,Ix_global);

function reconstruction = create_approx_solution_aux(W,given_sol,intermediate,N_proj,d_set,res,Ix_global)

x = intermediate;
n = length(x);
G = W;
Ix = Ix_global;
if size(d_set,1)<size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

% [x Ix] = round_numerically_eq(x,d_set,Ix); % it isnt necessary cause Ix_global is known
G(:,Ix) = []; % remove columns (unknowns)
Ix_c = setdiff([1:n],Ix); length(Ix_c)

rows = sum(G,2) < N_proj;
G(rows,:)=[]; % remove rows of G
rows0 = sum(G,2) == 0; % identify the null rows
G(rows0,:)=[]; % remove rows of G

% G=A

while ~isempty(G)
    
    X0 = rand(size(G,2),1);
    right_h_side = -G*X0;
    [Beta, ~, ~] = cgls_W(G, right_h_side, [], 1000, 1e-15);
    Beta = Beta+X0;
    res = res + norm(G*Beta,inf);
    
    if norm(Beta - zeros(length(Beta),1),1)<1e-14
        if norm(G,inf)>N_proj
            error('norm(G,inf)-N_proj=%d',norm(G,inf)-N_proj)
        end
        break
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
    
    if norm(W*(given_sol-x),inf) > (N_proj-1)*(d_set(end)-d_set(1))+res % MAKE RES MORE PRECISE
        error('3');
    end
    
    Ix = union(Ix,idx); % pixels already defined
    G = W;
    G(:,Ix) = []; % pixels removed from the linear system
    Ix_c = setdiff([1:n],Ix); length(Ix_c)
    
    rows = sum(G,2) < N_proj;
    G(rows,:)=[]; % remove rows of G
    rows0 = sum(G,2) == 0; % identify the null rows
    G(rows0,:)=[]; % remove rows of G
end
reconstruction = round2d_set(x, given_sol, d_set, Ix_c);

function [x Ix] = round_numerically_eq(x,d_set,Ix)
%
for i = setdiff([1:length(x)],Ix)
    for j = 1:length(d_set)
        if abs(x(i)-d_set(j)) <= 1e-15%1/sqrt(length(x))
            x(i) = d_set(j);
            Ix = union(Ix,i);
        end
    end
end