function x = cgls_pocs(A,b,x0,tol,d_set)

[x, res, sol] = cgls_W(A, b, x0, 10, 1e-5);
x = Pontod_set(x,d_set);
    
for iter = 20:10:90
    if norm(A*x-b,inf)<tol
        return
    end
    [x, res, sol] = cgls_W(A, b, x, iter, 1e-5);
    x = Pontod_set(x,d_set);
end
i = 9;
while norm(A*x-b,inf)>tol
    i = i+1;
    if mod(i,100)==0
        fprintf('%g\n',i)
    end
    [x, res, sol] = cgls_W(A, b, x, 100, 1e-5);
    x = Pontod_set(x,d_set);
end

function x = Pontod_set(x,d_set)

for i = 1:length(x)
    if x(i) > d_set(end)
        x(i) = d_set(end);
    elseif x(i) < d_set(1)
        x(i) = d_set(1);
    end
end