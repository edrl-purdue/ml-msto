% INTERPOLATE LOW VOLUME FRACTION CH
function CH = Interp_CH(ML_model, PTOMs, micro) %without material properties (E, k, etc)
    PTOMs0 = PTOMs;
    if micro.vf_cutoff > 0
        PTOMs(PTOMs0(:,1)<=micro.vf_cutoff,1) = micro.vf_cutoff;
    end
    CH = ML_model(PTOMs')';
    CH(PTOMs0(:,1)<=micro.vf_cutoff,:) = micro.Cminvec + (CH(PTOMs0(:,1)<=micro.vf_cutoff,:) - micro.Cminvec).*(PTOMs0(PTOMs0(:,1)<=micro.vf_cutoff,1)/micro.vf_cutoff).^1; % have not added this to piecewise ann
    % 3D Ortho Mecha: C11, C21, C22, C31, C32, C33, C44, C55, C66
    % to
    % 3D Aniso Mecha: C11, C21, C22, C31, C32, C33, C41, C42, C43, C44, C51, C52, C53, C54, C55, C61, C62, C63, C64, C65, C66
    if micro.dim == 3
        CH = [CH(:,1:6), zeros(size(CH,1), 3), CH(:,7), zeros(size(CH,1), 4), CH(:,8), zeros(size(CH,1), 5), CH(:,9)];
    end
end