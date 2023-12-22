% TOM TOPOLOGY OPTIMIZATION
function [xMicro, c_vec] = Micro_TO(scale,TOM_e,xMicro0,penal,maxloop,tolx,mat,KE)
    %% PREPARE MICROSCALE BOUNDARY CONDITIONS
    vf_e = TOM_e(1);
    weight_e = TOM_e(2);
    F0_Micro = cell(1,3); % KU = F for 3 physics
    U0_Micro = cell(1,3);
    N_Micro = cell(1,3);
    for p = 1:length(scale.phys)
        [F0_Micro{scale.phys(p)},U0_Micro{scale.phys(p)},N_Micro{scale.phys(p)}] = Prepare_BCMicro(scale,TOM_e,scale.phys(p));
    end
    
    %% PREPARE OC
    x = xMicro0;
    xMicro = x;
    change = 1;
    loop = 0;
    
    % SELECT OBJECTIVE FUNCTION
    % 1 = Compliance
    % 2 = Force-Displacement Error
    % 3 = Mutual Potential Energy
    micro_TO_obj = 1; % still need to sort out/organize the code for other objective functions
    
    %% MICROSCALE TOPOLOGY OPTIMIZATION LOOP
    while change > tolx && loop < maxloop
        loop = loop+1;
        
        %% MICROSCALE FEA
        [U_Micro,~,~] = FEA(xMicro,scale,penal,mat,KE,F0_Micro,U0_Micro,N_Micro);
        
        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        [c, dc, c_vec] = Obj_Fns(xMicro,scale,U_Micro,penal,mat,KE,[],[],weight_e); % Case 1: Calling Obj_Fns(xMicro,scale,U_Micro,penal,mat,KE,[],[],weight_e) is c and dc for SIMP stiffness matrix
        dv = ones(size(xMicro));
        dc = reshape(dc,size(xMicro));
        
        %% DENSITY FILTER
        dc(:) = scale.den_H*(xMicro(:).*dc(:))./scale.den_Hs./max(1e-3,xMicro(:));
        dc(:) = scale.den_H*(dc(:)./scale.den_Hs);
        dv(:) = scale.den_H*(dv(:)./scale.den_Hs);
        
        %% OPTIMALITY CRITERIA UPDATE
        l1 = 0; l2 = 1e9; move = 0.2;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            min_1 = min(x+move,x.*sqrt(-dc./dv/lmid));
            min_2 = min(1,min_1);
            max_1 = max(x-move,min_2);
            xnew = max(0,max_1);
            xMicro(:) = (scale.den_H*xnew(:))./scale.den_Hs;
            if sum(xMicro(:)) > vf_e*scale.nele
                l1 = lmid;
            else
                l2 = lmid;
            end
        end
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        
        %% PRINT RESULTS AND PLOT DENSITIES
        if scale.displayflag
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f TargetVol.:%7.3f ch.:%7.3f \n',loop,c,mean(xMicro(:)),vf_e,change);
            figure(3)
            clf;
            display_top(xMicro);
        end
    end
end