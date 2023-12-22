% ADD MESH INFORMATION TO A SCALE'S STRUCTURE
function scale = mesh_info(scale)
    scale.nele = max(scale.nelx*scale.nely,scale.nelx*scale.nely*scale.nelz); % number of elements
    if scale.dim == 2
        scale.nen = 4; % number of nodes per element
    elseif scale.dim == 3
        scale.nen = 8; % number of nodes per element
    end
end