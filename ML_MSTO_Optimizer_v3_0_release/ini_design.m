% PREPARE INITIAL DESIGN
function x0 = ini_design(scale,vf,lb,ub)
    beta = ub - lb; % range of distribution [0,1]

    if scale.design_ini == 1 % Uniformly Distributed Volume Fraction Design
        x0 = ones(scale.nelx,scale.nely,max(1,scale.nelz))*vf;
    elseif scale.design_ini == 2 % Normally Distributed Volume Fraction Design
        if vf < 0.5
            a = 0.0; b = 2*(vf - 0);
        elseif vf == 0.5
            a = 0.0; b = 1.0;
        elseif vf > 0.5
            a = 1 - 2*(1 - vf); b = 1.0;
        end
        range = ((b - a)/2)*(1 - beta); a = a + range; b = b - range;
        x0 = normrnd(vf,((b-a)/2)/4,scale.nelx,scale.nely,max(1,scale.nelz)); % Set to standard deviation of 4 (i.e., 99.994% of cells are within range)
        x0(x0<a) = a;
        x0(x0>b) = b;
    elseif scale.design_ini == 3 % Randomly Distributed Volume Fraction Design
        if vf < 0.5
            a = 0.0; b = 2*(vf - 0);
        elseif vf == 0.5
            a = 0.0; b = 1.0;
        elseif vf > 0.5
            a = 1 - 2*(1 - vf); b = 1.0;
        end
        range = ((b - a)/2)*(1 - beta); a = a + range; b = b - range;
        x0 = a + rand(scale.nelx,scale.nely,max(1,scale.nelz)).*(b-a);
    end
end