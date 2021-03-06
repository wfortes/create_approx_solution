function reconstruction = round2d_set(reconstruction, given_sol, d_set, Ix_c)
% round to available grey values
%
% Wagner Fortes 2014/2015 wfortes@gmail.com

d_set = sort(d_set);
for i = Ix_c
    if given_sol(i) >= d_set(end)
        reconstruction(i) = d_set(end);
    elseif given_sol(i) <= d_set(1)
        reconstruction(i) = d_set(1);
    else
        indice = find(d_set > given_sol(i),1); 
        if given_sol(i) <= (d_set(indice)+d_set(indice-1))/2;
            reconstruction(i) = d_set(indice-1);
        else
            reconstruction(i) = d_set(indice);
        end
    end
end