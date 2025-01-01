%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%          BARNAUD Rudy         %
%     Num Met 4 Phys - Ex 3.8   %
%           15 Oct 24           %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q] = Normalization(A)
    
    % Function normalizing vectors
    % Matrix A contains vectors (in row convention)
    
    Q = zeros(size(A));
    for k=1:size(A, 1)
        Q(k, :) = A(k, :)/norm(A(k, :));
    end
end