% CHECK CORRDINATES FOR POINT OF INTEREST
function output = isPoint(name,idnums,optional_ids,i,lcoorp,xyzc,dx,dy,dz)
% Checks if the current point (lcoorp(i), lcoorp(i+1), lcoorp(i+2)) is
% within ### percent of the target feature (name, idnums)
percent = 0.05; %                                                             NAME    FEATURE                                  ##
%                 DOMAIN NODE AND FACE NUMBERINGS                            'midp' = Midpoint of domain (point)/(area/volume) 00
%      2D Nodes               3D Nodes                  3D Faces             'corn' = Corner (point)/(area/volume)             2 x Corner #
%    3 ________ 4         5 ________ 7             ________                  'edge' = Edge (line)/(partial line)               Corner # to Corner #
%     |        |  ^+Y      |\       |\    ^+Z     |\   2   |\    1 = bottom  'face' = Face (surface)/(volume)                  2 x Face #
%     |  1 =   |  |        |6\ _____|_\ 8 |       | \ _____|_\   2 = top     'mide' = Midpoint of edge (point)/(line)          Corner # to Corner #
%     |  face  |  |        |  |     |  |  |       |3 | 5   |  |  3 = front   'midx' = x-Midline (line)/(area)                  2 x Face #
%    1|________|2 +---->  1|__|_____|3 |  +---->  |__|_____| 4|  4 = back    'midy' = y-Midline (line)/(area)                  2 x Face #
%                    +X     \ |      \ |   \  +Y   \ |   6  \ |  5 = left    'midz' = z-Midline (line)/(area)                  2 x Face #
%                           2\|_______\|4   \       \|_______\|  6 = right   'midf' = Midpoint of face (point)/(area)          2 x Face #
%                                            v +X        1                   'mids' = Surface btwn faces (surface)/(volume)    Face # to Face #
%                                                                                     (feature loaded/constrained)/(optional feature loaded/constrained)
%                                                                            OPTIONAL
%                                                                            '##_@@' where # is the feature number and @ is the percent of domain to cover. Usually ## is the same feature numbers specified in the second column. The exception to this is 'edge'
%                                                                            For example {'edge', 13, 12, '13_20'} as a mechanical support means that the middle 20% of the 13 edge should be fixed in the x and y directions
%                                                                                        {'edge', 13, 12, '11_20'} as a mechanical support means that the left (node 1 side) 20% of the 13 edge should be fixed in the x and y directions
    output = false;
    x_add_p = 0; y_add_p = 0; z_add_p = 0;
    if strcmp(name,'midp') % Midpoint of domain
        if idnums == 00
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
                if xyzc(1) == 1 % 3D
                    z_add_p = 0.5*str2double(optional_ids(4:end))/100;
                end
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'corn') % Corner
        if idnums == 11
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 22
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 33
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 44
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 55
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 66
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 77
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 88
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                y_add_p = str2double(optional_ids(4:end))/100;
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'edge') % Edge
        if idnums == 12 || idnums == 21
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 11
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 22
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                   (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 24 || idnums == 42
            if ~isempty(optional_ids)
                y_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 22
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 44
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                   (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 13 || idnums == 31
            if ~isempty(optional_ids)
                y_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 11
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 33
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                   (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 34 || idnums == 43
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 33
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 44
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                   (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 56 || idnums == 65
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 55
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 66
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                   (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 68 || idnums == 86
            if ~isempty(optional_ids)
                y_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 66
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 88
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                   (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 57 || idnums == 75
            if ~isempty(optional_ids)
                y_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 55
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 77
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                   (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 78 || idnums == 87
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 77
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 88
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                   (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 15 || idnums == 51
            if ~isempty(optional_ids)
                z_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 11
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 55
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                   (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymin
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 26 || idnums == 62
            if ~isempty(optional_ids)
                z_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 22
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 66
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymin
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                   (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymin
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 37 || idnums == 73
            if ~isempty(optional_ids)
                z_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 33
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 77
                    if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmin
                   (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymax
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 48 || idnums == 84
            if ~isempty(optional_ids)
                z_add_p = str2double(optional_ids(4:end))/100;
                if str2double(optional_ids(1:2)) == 44
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmin
                        output = true;
                    end
                elseif str2double(optional_ids(1:2)) == 88
                    if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                       (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymax
                       (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmax
                        output = true;
                    end
                end
            else
                if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmax
                   (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymax
                    output = true;
                end
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'face') % Face
        if idnums == 11
            if ~isempty(optional_ids)
                if xyzc(1) == 1 % 3D
                    z_add_p = str2double(optional_ids(4:end))/100;
                end
            end
            if (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p) %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 22
            if ~isempty(optional_ids)
                z_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p) %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 33
            if ~isempty(optional_ids)
                y_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) %ymin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 44
            if ~isempty(optional_ids)
                y_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) %ymax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 55
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p)   %xmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 66
            if ~isempty(optional_ids)
                x_add_p = str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p)   %xmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'mide') % Midpoint of edge
        if idnums == 12 || idnums == 21
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymin
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 24 || idnums == 42
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 13 || idnums == 31
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 34 || idnums == 43
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymax
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 56 || idnums == 65
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymin
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 68 || idnums == 86
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 57 || idnums == 75
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 78 || idnums == 87
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymax
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 15 || idnums == 51
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymin
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 26 || idnums == 62
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymin
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 37 || idnums == 73
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymax
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 48 || idnums == 84
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymax
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'midx') % x-Midline
        if idnums == 11
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 22
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 33
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymin
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 44
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymax
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 55 || idnums == 66
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'midy') % y-Midline
        if idnums == 11
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 22
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 33 || idnums == 44
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 55
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 66
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'midz') % z-Midline
        if idnums == 11 || idnums == 22
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 33
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p)                       %ymin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 44
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p)                       %ymax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 55
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 66
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p)       %ymid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'midf') % Midpoint of face
        if idnums == 11
            if ~isempty(optional_ids)
                if xyzc(1) == 1 % 3D
                    x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                    y_add_p = 0.5*str2double(optional_ids(4:end))/100;
                end
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(5)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmin
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 22
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - xyzc(6)) < (min([dx,dy,dz]*percent)) + z_add_p)                       %zmax
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 33
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(3)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymin
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 44
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p) &&...   %xmid
               (abs(lcoorp(i+1) - xyzc(4)) < (min([dx,dy,dz]*percent)) + y_add_p) &&...                 %ymax
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 55
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(1)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmin
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 66
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - xyzc(2)) < (min([dx,dy,dz]*percent)) + x_add_p) &&...                   %xmax
               (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) &&... %ymid
               (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p)       %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
    
    if strcmp(name,'mids') % Middle surface between two faces
        if idnums == 12 || idnums == 21
            if ~isempty(optional_ids)
                z_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+2) - mean([xyzc(5) xyzc(6)])) < (min([dx,dy,dz]*percent)) + z_add_p) %zmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 34 || idnums == 43
            if ~isempty(optional_ids)
                y_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i+1) - mean([xyzc(3) xyzc(4)])) < (min([dx,dy,dz]*percent)) + y_add_p) %ymid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        elseif idnums == 56 || idnums == 65
            if ~isempty(optional_ids)
                x_add_p = 0.5*str2double(optional_ids(4:end))/100;
            end
            if (abs(lcoorp(i) - mean([xyzc(1) xyzc(2)])) < (min([dx,dy,dz]*percent)) + x_add_p)   %xmid
                output = true;
            end
            x_add_p = 0; y_add_p = 0; z_add_p = 0;
        end
    end
end