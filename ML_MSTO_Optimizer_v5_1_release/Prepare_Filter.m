% PREPARE FILTER
function [H,Hs] = Prepare_Filter(scale,rmin)
    if scale.dim == 2 % 2D
        iH = ones(scale.nely*scale.nelx*((ceil(rmin)-1)+1)^2,1); jH = ones(size(iH)); sH = zeros(size(iH)); k = 0;
        for i1 = 1:scale.nely
            for j1 = 1:scale.nelx
                e1 = (i1-1)*scale.nelx+j1;
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),scale.nely)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),scale.nelx)
                        e2 = (i2-1)*scale.nelx+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                    end
                end
            end
        end
    elseif scale.dim == 3 % 3D
        iH = ones(scale.nelx*scale.nely*scale.nelz*(2*(ceil(rmin)-1)+1)^2,1); jH = ones(size(iH)); sH = zeros(size(iH)); k = 0;
        for k1 = 1:scale.nelz
            for i1 = 1:scale.nelx
                for j1 = 1:scale.nely
                    e1 = (k1-1)*scale.nelx*scale.nely + (i1-1)*scale.nely+j1;
                    for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),scale.nelz)
                        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),scale.nelx)
                            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),scale.nely)
                                e2 = (k2-1)*scale.nelx*scale.nely + (i2-1)*scale.nely+j2;
                                k = k+1;
                                iH(k) = e1;
                                jH(k) = e2;
                                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                            end
                        end
                    end
                end
            end
        end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);
end