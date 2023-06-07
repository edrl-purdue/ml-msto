% GENERALIZED OPTIMALITY CRITERION ALGORITHM
function [Xin, c] = alg_GOCM(FUN,Xin,A,B,Aeq,Beq,LB,UB,nonlcon,options,macro)
    % A more general OC implementation. Can handle objective functions and
    % sensitivities with positive or negative values. Can handle multiple
    % inequality constraints (A, B). Equality constraints (Aeq, Beq) are
    % converted to inequality constraints: Aeq*Xin <= Beq
    % nonlcon not considered
    
    % Prepare GOCM
    lmid = zeros(options.maxloop+1,size([Aeq; A],1));
    lmid(1,:) = [1, options.lmin*ones(1,size(A,1))]; % initial lagrange multiplier value (eq 11)
    A = [Aeq; A]; B = [Beq; B]; % converty equality constraints to inequality constraints
    dg = transpose(A./B); % sensitivity of linear, normalized inequality constraints
    dg = (macro.den_H*dg)./macro.den_Hs; % apply density filter
    change = 1;
    loop = 0;
    
    % Run GOCM
    while change > options.tolx && loop < options.maxloop
        loop = loop+1;

        % homogenized macroscale objective and derivative
        [c, dc] = FUN(Xin);
        if loop == 1; c0 = c; end % note initial compliance
        dc = dc/c0; % normalize sensitivies
    
        % GOC update
        g = transpose(A*Xin./B - 1); % normalized inequality constraints
        if loop == 1; g_change = zeros(size(g)); else; g_change = g-g_old; end % change in inequality constraints
        g_old = g;
        for i = 1:numel(g)
            if (g(i) > 0 && g_change(i) > 0) || (g(i) < 0 && g_change(i) < 0) % when the sign of g and g_change are the same then p0 is positive
                p0 = 1.0;
            elseif (g(i) > 0 && g_change(i) > -options.tolp) || (g(i) < 0 && g_change(i) < options.tolp)
                p0 = 0.5;
            else
                p0 = 0;
            end
            lmid(loop+1,i) = lmid(loop,i)*(1 + p0*(g(i) + g_change(i)));
        end
        lmid(loop+1,:) = min(max(lmid(loop+1,:),options.lmin),options.lmax);
        D = Xin.*(-(min(dc,0) + sum(min(lmid(loop+1,:).*dg,0),2))./(max(dc,0) + sum(max(lmid(loop+1,:).*dg,0),2))).^options.eta;
        min_1 = min(Xin+options.move, D); % minimum operation 1
        min_2 = min([UB, min_1], [], 2); % minimum operation 2
        max_1 = max(Xin-options.move, min_2); % maximum operation 1
        xnew = max([LB, max_1], [], 2); % maximum operation 2
        change = max(abs(xnew-Xin));
        Xin = xnew;
    end
end