function maskSet = blendedPolymask(curves, xVec, yVec, zVec)
%BLENDEDPOLYMASK   Creates a 3D mask "lofted" from a set of polygons
%   BW = BLENDEDPOLYMASK(C,X,Y,Z) returns a logical volume of size equal to
%   [length(Y) length(X) length(Z)]. C is a cell array of N-by-3 XYZ
%   coordinates specifying individual (planar) XY polygons. The first
%   Z-coordinate of each polygon in C specifies its planar location. X, Y
%   and Z are vectors specifying pixel locations in the 3D output volume.
%   The output BW is a "loft" from one polygon to the next. Interpolation
%   between polygons is done via the distance function of the masks at each
%   neighbouring polygon. The IP toolbox BWDIST function is required.
%
%   Example:
%     [circXY(:,1),circXY(:,2)] = pol2cart(linspace(0,2*pi,50)', 1);
%     sqXY = [-1 -1;1 -1;1 1;-1 1; -1 -1];
%     C = {[sqXY*5 ones(5,1)]           % Start with a small square
%         [circXY*40 ones(50,1)*30]     % Blend to a large circle
%         [sqXY*20 ones(5,1)*65]        % Blend to a large square
%         [circXY*10 ones(50,1)*99]};   % Blend to a small circle
%     X = linspace(-40, 40, 200);
%     Y = linspace(-40, 40, 200);
%     Z = linspace(0, 100, 400);
%     BW = blendedPolymask(C,X,Y,Z);
%     figure, patch(isosurface(X,Y,Z,BW,0.5),'FaceColor','g','EdgeColor','none','FaceAlpha',0.5)
%     view(3), camlight, hold on, axis image
%     cellfun(@(x)patch(x(:,1),x(:,2),x(:,3),'b'),C)
%
%   See also poly2mask

% Written by Sven Holcombe

% Initialise maskSet to the size of the input xVec, yVec, zVecs
lenX = length(xVec); lenY = length(yVec); lenZ = length(zVec);
maskSet = false([lenY lenX lenZ]);

% For each zVec location, exactly where along "curves" is it placed? (partial indices used)
numCurves = length(curves);
curveZlocations = cellfun(@(xyz)xyz(1,3),curves);
idxsIntoCurves = interp1(curveZlocations, 1:numCurves, zVec(:), 'linear');

% For any zSamps outside the range of curves, bring them back to the first/last curve
validSheetNos = find(~isnan(idxsIntoCurves));

% For each zVec, get its which are its nearest inferior/superior curve indices?
idxsLowHigh = [floor(idxsIntoCurves) ceil(idxsIntoCurves)];
[usedIdxs,~,validUsedGrpNo] = unique(idxsLowHigh(validSheetNos,:));

% Build BW masks and mask dist maps for each non-empty zLevel
usedBWmasks = arrayfun(@(i)poly2mask(...
    interp1(xVec, 1:lenX, curves{i}(:,1),'linear','extrap'),...
    interp1(yVec, 1:lenY, curves{i}(:,2),'linear','extrap'),...
    lenY, lenX) , usedIdxs, 'UniformOutput',false);

% For any curveZlocations exactly hitting a usedBWmask, use it directly
directCopySlices = find(diff(idxsLowHigh,[],2)==0); % Idxs into "usedIdxs"
for i = 1:length(directCopySlices)
    maskSet(:,:,directCopySlices(i)) = usedBWmasks{validUsedGrpNo(validSheetNos==directCopySlices(i))};
end

% Loop through remaining levels. If a pixel is closer to an "on"
% neighbour than an "off" neighbour, it gets turned on.
for i = setdiff(validSheetNos(:)', directCopySlices)
    twoDists = abs(idxsLowHigh(i,:) - idxsIntoCurves(i));
    twoInds = [find(usedIdxs==idxsLowHigh(i,1)) find(usedIdxs==idxsLowHigh(i,2))];
    BW1 = usedBWmasks{twoInds(1)}; BW2 = usedBWmasks{twoInds(2)};
    scaledDist1 = (bwdist(BW1) - bwdist(~BW1)) * twoDists(2);
    scaledDist2 = (bwdist(BW2) - bwdist(~BW2)) * twoDists(1);
    maskSet(:,:,i) = scaledDist1+scaledDist2 <= 0;
    % BELOW SHOWS COMPARISONS BETWEEN BW1, BW2, and interpolated result
    %         bb1 = bwboundaries(BW1,'noholes');  bb2 = bwboundaries(BW2,'noholes');  bb = bwboundaries(maskSet(:,:,i),'noholes');
    %         figure, plot(bb1{1}(:,2),bb1{1}(:,1),'r',bb2{1}(:,2),bb2{1}(:,1),'b',bb{1}(:,2),bb{1}(:,1),':g','LineWidth',2)
    %         title(sprintf('Rdist=%f, Bdist=%f',twoDists)), axis image; uiwait
end
