% HR 12/03/22 PDA function
function [D,w,m] = product_difference_algorithm(x,PDF,N)
    for k = 1:(2*N)
        m(k) = sum(PDF(:).*x(:).^(k-1));
    end
    
    P = zeros(2*N+1);
    P(1,1) = 1;
    for i = 1:(2*N)
        P(i,2) = (-1)^(i-1) * m(i);
    end
    
    for j = 3:(2*N+1)
        for i = 1:(2*N+2-j)
            P(i,j) = P(1,j-1)*P(i+1,j-2) - P(1,j-2)*P(i+1,j-1);
        end
    end

    alpha(1) = 0;
    for i = 2:2*N
        alpha(i) = P(1,i+1)/(P(1,i)*P(1,i-1));
    end
    
    for i = 1:N
        a(i) = alpha(2*i) + alpha(2*i - 1);
    end
    
    for i = 1:(N-1)
        b(i) = -sqrt(alpha(2*i + 1)*alpha(2*i));
    end
    
    J = gallery('tridiag',b,a,b)
    [w,D] = eig(full(J))
    
    D = diag(D);
%     w = w(1,:).^2;
%     w = w';
    w = w(:,1).^2

end