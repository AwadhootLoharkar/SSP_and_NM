%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%     Num Met 4 Phys - Ex 3.8   %
%           15 Oct 24           %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R] = GramSchmidtClassic(A)
    
    % Gram-Schmidt process as it's done usually
    % This is not a stable process numerically
    % Matrix A contains vectors (in row convention)
    
    if(size(A, 1) ~= size(A, 2))
        error("Input matrix is not square.")
    end
    
    %% Equation 3.120 & 3.119
    Q = zeros(size(A));
    for k=1:size(A, 1)
        Q(k, :) = A(k, :);
        for j=1:k-1
            Q(k, :) = Q(k, :) - dot(Q(j, :), A(k, :))/dot(Q(j, :), Q(j, :))*Q(j, :);
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