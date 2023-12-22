% CALCULATE SHAPE FUNCTION AND ITS DERIVATIVE
% ncoord = [xi, eta, zeta] vector of local/natural coordinates
function [N, dN] = N_dN_calc(ncoord)
    if length(ncoord) == 2
        xi = ncoord(1); eta = ncoord(2);
        % Shape function
        N(1) = 0.25*(1.0 - xi)*(1.0 - eta);
        N(2) = 0.25*(1.0 + xi)*(1.0 - eta);
        N(3) = 0.25*(1.0 + xi)*(1.0 + eta);
        N(4) = 0.25*(1.0 - xi)*(1.0 + eta);
        %  With respect to xi:
        dN(1,1) = -0.25 * (1.0 - eta);
        dN(2,1) = 0.25 * (1.0 - eta);
        dN(3,1) = 0.25 * (1.0 + eta);
        dN(4,1) = -0.25 * (1.0 + eta);
        %  With respect to eta:
        dN(1,2) = -0.25 * (1.0 - xi);
        dN(2,2) = -0.25 * (1.0 + xi);
        dN(3,2) = 0.25 * (1.0 + xi);
        dN(4,2) = 0.25 * (1.0 - xi);
    elseif length(ncoord) == 3
        xi = ncoord(1); eta = ncoord(2); zeta = ncoord(3);
        % Shape function
        N(1) = 0.125*(1.0 - xi)*(1.0 - eta)*(1.0 - zeta);
        N(2) = 0.125*(1.0 + xi)*(1.0 - eta)*(1.0 - zeta);
        N(3) = 0.125*(1.0 + xi)*(1.0 + eta)*(1.0 - zeta);
        N(4) = 0.125*(1.0 - xi)*(1.0 + eta)*(1.0 - zeta);
        N(5) = 0.125*(1.0 - xi)*(1.0 - eta)*(1.0 + zeta);
        N(6) = 0.125*(1.0 + xi)*(1.0 - eta)*(1.0 + zeta);
        N(7) = 0.125*(1.0 + xi)*(1.0 + eta)*(1.0 + zeta);
        N(8) = 0.125*(1.0 - xi)*(1.0 + eta)*(1.0 + zeta);
        %  With respect to xi:
        dN(1,1) = -0.125 * (1.0 - eta) * (1.0 - zeta);
        dN(2,1) = 0.125 * (1.0 - eta) * (1.0 - zeta);
        dN(3,1) = 0.125 * (1.0 + eta) * (1.0 - zeta);
        dN(4,1) = -0.125 * (1.0 + eta) * (1.0 - zeta);
        dN(5,1) = -0.125 * (1.0 - eta) * (1.0 + zeta);
        dN(6,1) = 0.125 * (1.0 - eta) * (1.0 + zeta);
        dN(7,1) = 0.125 * (1.0 + eta) * (1.0 + zeta);
        dN(8,1) = -0.125 * (1.0 + eta) * (1.0 + zeta);
        %  With respect to eta:
        dN(1,2) = -0.125 * (1.0 - xi) * (1.0 - zeta);
        dN(2,2) = -0.125 * (1.0 + xi) * (1.0 - zeta);
        dN(3,2) = 0.125 * (1.0 + xi) * (1.0 - zeta);
        dN(4,2) = 0.125 * (1.0 - xi) * (1.0 - zeta);
        dN(5,2) = -0.125 * (1.0 - xi) * (1.0 + zeta);
        dN(6,2) = -0.125 * (1.0 + xi) * (1.0 + zeta);
        dN(7,2) = 0.125 * (1.0 + xi) * (1.0 + zeta);
        dN(8,2) = 0.125 * (1.0 - xi) * (1.0 + zeta);
        %  With respect to zeta:
        dN(1,3) = -0.125 * (1.0 - xi) * (1.0 - eta);
        dN(2,3) = -0.125 * (1.0 + xi) * (1.0 - eta);
        dN(3,3) = -0.125 * (1.0 + xi) * (1.0 + eta);
        dN(4,3) = -0.125 * (1.0 - xi) * (1.0 + eta);
        dN(5,3) = 0.125 * (1.0 - xi) * (1.0 - eta);
        dN(6,3) = 0.125 * (1.0 + xi) * (1.0 - eta);
        dN(7,3) = 0.125 * (1.0 + xi) * (1.0 + eta);
        dN(8,3) = 0.125 * (1.0 - xi) * (1.0 + eta);
    end
end