% INTERPOLATE LOW VOLUME FRACTION CH
function CH = Interp_CH(ML_model, PTOMs, micro) %without material properties (E, k, etc)
    PTOMs0 = PTOMs;
    if micro.vf_cutoff > 0
        PTOMs(PTOMs0(:,1)<=micro.vf_cutoff,1) = micro.vf_cutoff;
    end
    CH = ML_model(PTOMs')';
    CH(PTOMs0(:,1)<=micro.vf_cutoff,:) = micro.Cminvec + (CH(PTOMs0(:,1)<=micro.vf_cutoff,:) - micro.Cminvec).*(PTOMs0(PTOMs0(:,1)<=micro.vf_cutoff,1)/micro.vf_cutoff).^1; % have not added this to piecewise ann
end