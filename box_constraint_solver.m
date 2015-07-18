function x = box_constraint_solver(A,b,maxit,tol,method)

x = zeros(size(A,2),1);
x_aux = x;
iter = 0;

if strcmp(method,'ART')
%     [x_ls, res, sol] = cgls_W(A, b, [], 100, 1e-10);
%     x = Pontod_set(x_aux,[0;1]);
    while norm(b-A*x,inf) > tol && iter < maxit
        for i = 1:length(b)
            if norm(A(i,:))>1e-15
                c = (b(i)-dot(A(i,:),x))/norm(A(i,:))^2;
                x_aux = x_aux+c*A(i,:)';
                if mod(iter,100)==0
                    fprintf('%d\t%g\n',iter,norm(b-A*x,inf))
                end
                x = Pontod_set(x_aux,[0;1]);
                iter = iter+1;
            end
        end
    end
elseif strcmp(method,'ART-tol')
    while norm(b-A*x,inf) > tol && iter < maxit
        for i = 1:length(b)
            if norm(A(i,:))>1e-15
                c = (b(i)-dot(A(i,:),x))/norm(A(i,:))^2;
                x_aux = x_aux+c*A(i,:)';
                if mod(iter,100)==0
                    fprintf('%d\t%g\n',iter,norm(b-A*x,inf))
                end
                x = Pontod_set(x_aux,[0;1]);
                iter = iter+1;
            end
        end
    end
elseif strcmp(method,'SIRT')
    C = diag(sum(A,1));
    R = diag(sum(A,2));
    M = C*A'*R;
    while norm(b-A*x,inf) > tol && iter < maxit
        r = b - A*x;
        x_aux = x_aux+M*r;
        x = Pontod_set(x_aux,[0;1]);
        iter = iter+1;
    end
end

function x = Pontod_set(x,d_set)

for i = 1:length(x)
    if x(i) > d_set(end)
        x(i) = d_set(end);
    elseif x(i) < d_set(1)
        x(i) = d_set(1);
    end
end