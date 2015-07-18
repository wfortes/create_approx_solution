function reconstruction = create_approx_solution_piecewise3(W,given_sol,N_proj,d_set,key,thresh)

big_piece = sqrt(length(given_sol));
if strcmp(key,'divide')
    reconstruction = intermediate_solution(W,given_sol,given_sol,N_proj,d_set,big_piece/4,thresh);
%     reconstruction = intermediate_solution(W,given_sol,reconstruction,N_proj,d_set,big_piece/4,thresh);
    reconstruction = intermediate_solution(W,given_sol,reconstruction,N_proj,d_set,big_piece/2,thresh);
    reconstruction = intermediate_solution(W,given_sol,reconstruction,N_proj,d_set,big_piece,thresh);
else
    reconstruction = intermediate_solution(W,given_sol,given_sol,N_proj,d_set,big_piece,thresh);
end

function reconstruction_partial = intermediate_solution(W,given_sol,x_previous,N_proj,d_set,piece,thresh)

% this is only implemented for square images and divisible by the piece

n = length(given_sol);
mm = n/piece^2;
m = sqrt(mm);
reconstruction_partial = x_previous;
Ix_global = [];
res = [];

if size(d_set,1)<size(d_set,2);
    d_set = d_set';
end
d_set = sort(d_set);

[x Ix] = round_numerically_eq(x_previous,d_set,[],thresh);% global
Ix_global = union(Ix_global,Ix);

for i_pieces_columns = 1:m
    for i_pieces_rows = 1:m
        Ix = [];
        pixels_in_piece = [];
        
        for i_row = 0:piece-1
            pixels_in_piece_aux = i_row*piece*m+1+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1):i_row*piece*m+piece+(i_pieces_rows-1)*piece+piece^2*m*(i_pieces_columns-1);
            pixels_in_piece = union(pixels_in_piece,pixels_in_piece_aux); %global coordinates
        end
        x = x_previous(pixels_in_piece);
%         if piece^2 == n && norm(x-x_previous,2)~=0
%             error('x~=given_sol')
%         end
        
        [x Ix] = round_numerically_eq(x,d_set,Ix,thresh); % local
        Ix_c = setdiff(1:piece^2,Ix);
        
        G = W;
        G(:,Ix_global) = []; % remove columns from W
        rows = sum(G,2) < N_proj;
        
        A = W(:,pixels_in_piece); % remove columns from W
        A(:,Ix) = []; % remove columns local
        
        A(rows,:)=[]; % remove rows of A according to the set G
        rows0 = sum(A,2) == 0; % identify the null rows
        A(rows0,:)=[]; % remove null rows of A
        
        while ~isempty(A)
            
            X0 = rand(size(A,2),1);
            right_h_side = -A*X0;
            [Beta, ~, ~] = cgls_W(A, right_h_side, [], 1000, 1e-15);
            Beta = Beta+X0;
            res = res + norm(A*Beta,inf);
            if norm(Beta - zeros(length(Beta),1),1)<1e-14
%                 warning('null solution')
                break
            end
            
            Beta_aux = zeros(piece^2,1);
            Beta_aux(Ix_c) = Beta; % local
            Beta = Beta_aux;
            M_beta = diag(1./Beta);
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
            
            if isempty(thresh)
                reconstruction_partial(pixels_in_piece) = x;
            else
                [x Ix] = round_numerically_eq(x,d_set,Ix,thresh);
                reconstruction_partial(pixels_in_piece) = x;
            end
            
%             H=W;
%             H(rows,:)=[];
%             for au = 1:size(H,1)
%                 if abs(H(au,:)*(reconstruction_partial-given_sol))>res % MAKE RES MORE PRECISE
%                     error('1');
%                 end
%             end
%             HH=W;
%             HH(~rows,:)=[];
%             for au = 1:size(HH,1)
%                 if abs(HH(au,:)*(reconstruction_partial-given_sol))> (N_proj-1)*(d_set(end)-d_set(1))+res % MAKE RES MORE PRECISE
%                     error('2');
%                 end
%             end
            
            Ix = union(Ix,idx); % pixels already defined
            Ix_c = setdiff(1:piece^2,Ix); 
            %length(Ix_c)
            
            Ix_global = union(Ix_global,pixels_in_piece(idx));
            G = W;
            G(:,Ix_global) = []; % pixels removed from the set G
            rows = sum(G,2) < N_proj;
            
            A = W(:,pixels_in_piece); % remove columns of A
            A(:,Ix) = []; % remove columns of A (local)
            
%             if piece^2 == n && norm(A-G,2)~=0
%                 error('A~=G')
%             end
            
            A(rows,:)=[]; % remove rows of A according to the set G
            rows0 = sum(A,2) == 0; % identify the null rows
            A(rows0,:)=[]; % remove null rows of A
        end
        reconstruction_partial(pixels_in_piece) = x;
    end
end

if piece^2 == n
    reconstruction_partial = round2d_set(reconstruction_partial, given_sol, d_set, Ix_c);
end

function [x Ix] = round_numerically_eq(x,d_set,Ix,thresh)
%
if isempty(thresh)
    thresh = 1e-14;
end
for i = setdiff([1:length(x)],Ix)
    for j = 1:length(d_set)
        if abs(x(i)-d_set(j)) <= thresh
            x(i) = d_set(j);
            Ix = union(Ix,i);
        end
    end
end