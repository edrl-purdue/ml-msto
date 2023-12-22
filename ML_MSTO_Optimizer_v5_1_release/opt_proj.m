function abs_error = opt_proj(scaler, vf, scale, Mesh_line_data)
    Mesh_line_data(:,9:10) = Mesh_line_data(:,9:10)*scaler;
    Mesh_line_data(Mesh_line_data(:,9) < scale.thick_tol,9) = 0; Mesh_line_data(Mesh_line_data(:,10) < scale.thick_tol,10) = 0; % identify thin lattices
    Mesh_line_data(any(Mesh_line_data(:,9:10) < scale.thick_tol,2),:) = []; % remove thin lattices
    [Mesh_uni,~,ii] = unique(Mesh_line_data(:,1:2)); % identify disconnected lattices
    Mesh_line_data(any(Mesh_line_data(:,1) == Mesh_uni(accumarray(ii,1).' == 1)',2),:) = []; % remove disconnected lattices
    Mesh_line_data(any(Mesh_line_data(:,2) == Mesh_uni(accumarray(ii,1).' == 1)',2),:) = []; % remove disconnected lattices
    xProj = Project_SDF(scale, Mesh_line_data);
    abs_error = abs(mean(xProj(:)) - vf);
end