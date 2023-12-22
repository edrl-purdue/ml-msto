% MACROSCALE DENSITY GRADING CONSTRAINT
function [A, b] = density_grad(macro)
    if macro.delta == 1
        A = [];
        b = [];
    else
        Dmax = 1/macro.delta; % constraint range is greater than 1
        Dmax = 1; % constraint range is only immediately adjacent elemnts (e.g., NSEW in 2D)
        xNums = reshape(1:macro.nele,[macro.nelx, macro.nely, max(macro.nelz,1)]); % element nums at point a
        
        theta1 = []; theta2 = [];
        for e = 1:macro.nele
            dists = bwdist(xNums==e);
            theta1 = [theta1; e*ones(size(xNums(dists<=Dmax)))];
            theta2 = [theta2; xNums(dists<=Dmax)];
        end
        A = sparse(transpose(1:numel(theta1(:))),    theta1(:),     ones(size(theta1(:))),    numel(theta1(:)), macro.nele) +... % Theta_j
            sparse(transpose(1:numel(theta2(:))),    theta2(:),    -ones(size(theta2(:))),    numel(theta2(:)), macro.nele); % -Theta_k

        A = A(any(A,2),:); % Theta_j - Theta_k
        b = ones(size(A,1),1)*macro.delta;
    end
end