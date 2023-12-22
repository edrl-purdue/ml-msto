% OPTIMALITY CRITERION ALGORITHM
function [Xin, c] = alg_OCM(FUN,Xin,A,B,Aeq,Beq,LB,UB,nonlcon,options)
    % Not a general OC implementation. Only valid for a positive objective
    % function, always negative sensitivities, and one equality constraint.
    % The inequality constraints (A, B) are not considered
    % nonlcon not considered

    % Prepare OCM
    change = 1;
    loop = 0;

    % Run OCM
    while change > options.tolx && loop < options.maxloop
        loop = loop+1;

        % homogenized macroscale objective and derivative
        [c, dc] = FUN(Xin);
    
        % OC update
        l1 = options.lmin; l2 = options.lmax;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            D = Xin.*(-dc/lmid).^options.eta; 
            min_1 = min(Xin+options.move, D); % minimum operation 1
            min_2 = min([UB, min_1], [], 2); % minimum operation 2
            max_1 = max(Xin-options.move, min_2); % maximum operation 1
            xnew = max([LB, max_1], [], 2); % maximum operation 2
            if 0 > Beq - Aeq*xnew
                l1 = lmid;
            else
                l2 = lmid;
            end
        end
        change = max(abs(xnew-Xin));
        Xin = xnew;
    end
end