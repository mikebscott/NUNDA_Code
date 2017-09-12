function [cutCoords,vertlimit] = ascendingAortaLocation(labeledskel,mask_struct,varargin)
% Description

if nargin == 2
    plotflag = false;
elseif nargin == 3
    path = varargin{1};
    plotflag = true;
else
    error('Too many inputs')
end

% Grab the aorta mask
aortaMask = mask_struct.mask;
% Find the voxel size
vox = mask_struct.vox;
% Flatten the mask
flatMask = any(aortaMask,3);

% Flatten the labeled skel
flatSkel = max(labeledskel,[],3);
% Find the indices of the aorta midline
midlineIdx = find(flatSkel == 1);
% Convert to coordinates
[Xmidline,Ymidline] = ind2sub(size(flatSkel),midlineIdx);
% Find the indicies of the brachiocephalic artery
bcaIdx = find(flatSkel == 2);
% Convert to coordinates
[Xbca,Ybca] = ind2sub(size(flatSkel),bcaIdx);

% Find the X coordinate of the lowest value of the brachiocephalic branch
junctionX = max(Xbca)+1;
% Find the range of coordinates that are within the aorta at this level
idx = find(flatMask(junctionX,:) == 1);
minY = min(idx);
maxY = max(idx);
% Preallocate a vector to store the results
ycoords = minY:1:maxY;
distance = NaN(length(ycoords),1);
% Find the distance above the junction (perpendicular) that each aorta
% point is at
for ii = 1:length(ycoords)
    for jj = 0:junctionX-1
       if(flatMask(junctionX-jj,ycoords(ii)) == 0)
          distance(ii) = jj; 
          break;
       end
    end
end

% Find the difference in each distance point
distanceDiff = diff(distance);

% Find the coordinate of the first big jump (>=3 vox)
firstfallidx = find(distanceDiff < -2,1);
if(isempty(firstfallidx))
    firstfallidx = round(length(distanceDiff)/3);
end
if(max(distanceDiff > 3) && find(distanceDiff > 3,1) < firstfallidx)
    jumpIdx = find(distanceDiff > 3,1);
else
    % Testing: use the mean value close before firstfallidx
    testval = mean(distance(1:firstfallidx));
    testpoints = distance(1:firstfallidx);
    [~,jumpIdx] = min(abs(testpoints - testval));
    
    
%     possibleJumpIdx = find(distanceDiff == max(distanceDiff));
%     if(length(possibleJumpIdx) > 1)
%         if(length(possibleJumpIdx) == 2 && (possibleJumpIdx(1) == 1 || possibleJumpIdx(1) == 1))    
%                 jumpIdx = possibleJumpIdx(2);
%         else
%             
%             % Find the jumpIdx closest to 1/4 of the way across distanceDiff
%             goalIdx = round(length(distanceDiff)/4);
%             [~,temp] = min(abs(possibleJumpIdx - goalIdx));
%             jumpIdx = possibleJumpIdx(temp);
%             clear temp;
%         end
%     else
%         [~,jumpIdx] = max(distanceDiff);
%     end
end
Xcut = junctionX-distance(jumpIdx)+1;
Ycut = minY+jumpIdx-1;



% Find the X,Y,Z coordinates of the cut
Zcut = find(aortaMask(Xcut,Ycut,:) == 1);
Zcut = round(mean(Zcut));
% Find the coordinates of the midline
midlineIdx = find(labeledskel == 1);
[Xmidline,Ymidline,Zmidline] = ind2sub(size(aortaMask),midlineIdx);
% Convert to mm
XmidlineC = Xmidline .* vox(1);
YmidlineC = Ymidline .* vox(2);
ZmidlineC = Zmidline .* vox(3);
% Convert the cut coordinates to mm
Xcutmm = Xcut * vox(1);
Ycutmm = Ycut * vox(2);
Zcutmm = Zcut * vox(3);
% Calculate the distance from each midline point to the cut coordinate
distance = ((XmidlineC-Xcutmm).^2 + (YmidlineC-Ycutmm).^2 + (ZmidlineC-Zcutmm).^2).^0.5;
% Find the idx of the closest midline point
[~,closestIdx] = min(distance);

% Set the output
cutCoords = [XmidlineC(closestIdx) YmidlineC(closestIdx) ZmidlineC(closestIdx)];
vertlimit = Xcut;
% Plot and save
if(plotflag == true)
    fig101 = figure(101);
    imshow(flatMask)
    hold on;
    plot(Ymidline,Xmidline,'.r')
    plot(Ybca,Xbca,'.b')
    plot([Ycut Ycut],[1 size(flatMask,1)],'-g')
    plot([1 size(flatMask,2)], [Xcut Xcut], '-g')
    plot(Ymidline(closestIdx),Xmidline(closestIdx),'ok')
    saveas(fig101,[path filesep() 'AAcut.png'])
end
end