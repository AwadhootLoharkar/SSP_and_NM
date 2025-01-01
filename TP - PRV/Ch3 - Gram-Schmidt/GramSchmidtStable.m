%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%     Num Met 4 Phys - Ex 3.8   %
%           15 Oct 24           %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R] = GramSchmidtStable(A)
    
    % Stable Gram-Schmidt process
    % Matrix A contains vectors (in row convention)
    
    if(size(A, 1) ~= size(A, 2))
        error("Input matrix is not square.")
    end
    
    %% Equation 3.128 & 3.119
    Q = A;
    for k=1:size(A, 1)-1
        for i=(k+1):size(A, 1)
            Q(i, :) = Q(i, :) - dot(Q(k, :), Q(i, :))/dot(Q(k, :), Q(k, :))*Q(k, :);
        end
    end

    %% Normalization
    Q = Normalization(Q);

    %% Equation 3.135
    R = zeros(size(A));
    for i=1:size(A, 1)
        for j=i:size(A, 1)
            R(i, j) = dot(Q(i, :), A(j, :));
        end
    end
end