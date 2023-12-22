% NORMALIZE TOM BOUNDARY CONDITION VALUES
function Xe_int_Macro = Normalize_TOM_BC(Xe_int_Macro,scale)
    % Xe_int_Macro is an input vector of unnormalized macroscale elemental nodal field variables for the field variable X (displacements, forces, etc). These values are supplied as boundary conditions inputs to a subproblem
    
    for e = 1:size(Xe_int_Macro,1)
        ndofn = size(Xe_int_Macro,2)/scale.nen;
        if ndofn > scale.dim
            ndofn = scale.dim;
        end
        ssum = 0;
        for d = 1:ndofn
            ssum = ssum + Xe_int_Macro(e,d:ndofn:end).^2;
        end
        Xmag = max(ssum.^0.5);
        Xe_int_Macro(e,:) = Xe_int_Macro(e,:)/Xmag;
    end
end