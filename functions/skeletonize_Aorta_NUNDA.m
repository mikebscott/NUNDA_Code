function [labeledskel, interpolatedMidline, nodefinal, linkfinal] = skeletonize_Aorta_NUNDA(mrstruct_mask)
% skeletonize_Aorta.m
%
%   This code takes a mask of the aorta and finds the centerline and the
%   three major branches if possible. If this is not possible (due to loops
%   in the skeleton), the code will find the aorta and what it thinks is
%   the brachiocephalic artery.
%
%   Inputs: see example usage, below....
%   - mrstruct_mask: a structure that contains, at a minimum:
%       - mrstruct_mask.mask: a binary array that serves as a mask for the
%       region of interest
%       - mrstruct_mask.vox: a vector of length at least 3 that gives the
%       voxel dimensions in x,y,z (later dimensions are not used)
%
%   Constants:
%   THR: threshold for branch size, all smaller branches are pruned.
%   Default value is 3. Code runs faster with higher values usually, but
%   may miss branches. Can increase if you have nice geometry and good
%   segmentations!
%
%   Outputs:
%   - labeledskel: an x*y*z matrix with zero at every point, except the
%   skeleton. 1 indicates aorta, 2 brachiocephalic artery, 3 left common
%   carotid artery, 4 left subclavian artery
%   - interpolatedMidline: an nx3 array containing the coordinates of the
%   midline, in mm
%   - nodefinal: a structure containing information about the nodes.
%       .idx: linear index of the voxel
%       .links: links attached to the node
%       .conn: nodes that are connected by one link
%       .comx: x coordinate (voxel, not mm) of the node
%       .comy: y coordinate (voxel, not mm) of the node
%       .comz: z coordinate (voxel, not mm) of the node
%       .ep: 1 if the node is an end-point, 0 if not
%       .vessel: 1 = aorta, 2 = brachiocephalic, 3 = left common carotid,
%           4 = left subclavian
%   - linkfinal: a structure containing information about the links
%       .n1: number of the first node the link connects to
%       .n2: number of the second node the link connects to
%       .point: linear indicies of points in the link
%       .vessel: 1 = aorta, 2 = brachiocephalic, 3 = left common carotid,
%           4 = left subclavian
%
%   Example usage:
%   [labeledskel, interpolatedMidline, nodefinal, linkfinal] = skeletonize_Aorta_NUNDA(mrstruct_mask)
%
%   Required functions:
%   - skel2graph3d: from the Mathworks file exchange. Generates a graph
%   from the skeleton data
%   - Skeleton3D: from the Mathworks file exchange. Generates a skeleton of
%   a 3D volume. Implemented in Matlab by Philip Kollmannsberger, but see
%   the download page for references.
%   - mimics_to_mrstruct: located in the matlab_nu folder, converts the
%   mimics text files into an mrstruct that is easier to work with in
%   Matlab.
%
%   References:
%
%   Known bugs/shortcomings:
%   - Could add code to extract the left subclavian artery in addition to
%   the brachiocephalic in cases where there are not exactly three branches
%   found. Code would be identical to the code to extract the
%   brachiocephalic except switching min to max for the coordinates. Would
%   probably also be good to check that there are at least two
%   possibleBranchEnds and possibleBranchNodes as well, and that the same
%   branch hasn't been labeled as both the brachiocephalic and left
%   subclavian.
%   - Currently the outputs are not maximally simplified, they retain nodes
%   that used to have branches and links between them.
%
% Written by Mike Scott, August 2017 (Northwestern University)
% michael.scott1@northwestern.edu

%% Input check
% Default do not show plots
plotflag = false;

% Convert mask to binary
if isfield(mrstruct_mask, 'dataAy')
    aortaMask = (mrstruct_mask.dataAy ~= 0);
elseif isfield(mrstruct_mask, 'mask')
    aortaMask = (mrstruct_mask.mask ~= 0);
else
    error('Expecting the mask struct to have either a dataAy or mask field.')
end

% Define voxel sizes
voxelSize = mrstruct_mask.vox;

%% Process the skeleton
% Run skeleton3d
skel = Skeleton3D(aortaMask);
% Test
THR = 3;
[A,node,link] = Skel2Graph3D(skel,THR);

% Eliminate unconnected nodes
w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

numNodes = size(node,2);
newNumNodes = 0;
while(newNumNodes < numNodes)
    numNodes = size(node,2);
    skel = Graph2Skel3D(node,link,w,l,h);
    [A,node,link] = Skel2Graph3D(skel,THR);
    newNumNodes = size(node,2);
end

%% Find the longest path
% Preallocate a structure to store the following for each endpoint node:
%   Longest path (list of nodes)
%   Longest path length
% For each endpoint node
longestPathDistance = 0;
longestPathVector = [];
for ii = 1:size(node,2)
    if(node(ii).ep == 1)
        % Allocate a structure to hold the search data
        searchData = struct('path',{},'distance',{});
        % Call a recursive function to look at the connected links
        [searchData] = calculatePaths(node, searchData, ii);
        
        % Calculate the distance of each of the paths
        searchData = calculatePathlengths(searchData,full(A));
        
        % Extract only the longest path
        [longdist,longIdx] = max(cat(1,searchData.distance));
        % Compare to previous longest path
        if(longdist > longestPathDistance)
            % Save the distance and path vector
            longestPathDistance = searchData(longIdx).distance;
            longestPathVector = searchData(longIdx).path;
        end
        % Erase the search data structure
        clear searchData;
    end
end

%% Order the longest path from the aortic root to the descending aorta
if(node(longestPathVector(1)).comx > node(longestPathVector(end)).comx)
    % Reverse the order of the path if the first node is lower than the last
    longestPathVector = fliplr(longestPathVector);
end

%% Ensure the descending aorta is followed correctly
% To do this: find the ultimate branch point, then check which of the two
% branches has a smaller angle between it and the proximal section of the
% descending aorta.
% Check if the node with the highest y value is the one in the longest path
if(node(longestPathVector(end)).comx ~= max(cat(1,node.comx)))
    % Find the node with the highest x coordinate
    [~,testNode] = max(cat(1,node.comx));
    
    % First check if the testNode is already a part of longestPathVector
    if(any(longestPathVector == testNode))
        idx = find(longestPathVector == testNode);
        longestPathVector = longestPathVector(1:idx);
        clear idx;
    elseif isempty(node(testNode).links)
        % Do nothing, the testNode is not connected to the network
    else
        
        % Find the common branch point (this is more complicated than
        % anticipated, since both of the branches could have other branches
        % before sharing a common node.
        % Start by finding the shortest path between the nodes
        P = shortestpath(graph(A),testNode,longestPathVector(end));
        
        % Find the common branch point
        commonNode = 0;
        for ii = 2:length(P)
            if(any(longestPathVector == P(ii)))
                commonNode = P(ii);
                break;
            end
        end
        
        % Get the nodes closest to the branch point in the original longest and
        % test branches
        originalNextNode = longestPathVector(find(longestPathVector == commonNode)+1);
        testNextNode = P(find(P == commonNode)-1);
        
        % Get the points of the current longest branch, then fit a line to them
        %for ii = 1:3
        for ii = 1:length(node(commonNode).links)
            if(link(node(commonNode).links(ii)).n1 == commonNode)
                % If the common node is listed first (ie the common node is the
                % source of the link)
                % Get the first five points of the link
                if(length(link(node(commonNode).links(ii)).point) >= 5)
                    firstFive = link(node(commonNode).links(ii)).point(1:5);
                else
                    firstFive = link(node(commonNode).links(ii)).point;
                end
            else
                % If the common node is listed second (ie the common node is the
                % sink of the link)
                if(length(link(node(commonNode).links(ii)).point) >= 5)
                    firstFive = fliplr(link(node(commonNode).links(ii)).point(end-4:end));
                else
                    firstFive = fliplr(link(node(commonNode).links(ii)).point);
                end
            end
            
            % Check if the other node is the proximal descending aorta,
            % current longest branch, or branch being tested
            if(node(commonNode).conn(ii) == testNextNode)
                % If the test node
                % Convert the first and last points to coordinates
                [X,Y,Z] = ind2sub(size(aortaMask),firstFive);
                xu = X(end)-X(1);
                yu = Y(end)-Y(1);
                zu = Z(end)-Z(1);
                testNodeUnitVector = 1/sqrt(xu^2 + yu^2 + zu^2) * [xu yu zu];
            elseif(node(commonNode).conn(ii) == originalNextNode)
                % If the current longest branch node
                % Convert the first and last points to coordinates
                [X,Y,Z] = ind2sub(size(aortaMask),firstFive);
                xu = X(end)-X(1);
                yu = Y(end)-Y(1);
                zu = Z(end)-Z(1);
                longestPathUnitVector = 1/sqrt(xu^2 + yu^2 + zu^2) * [xu yu zu];
            else
                % If the points are in the proximal descending aorta, need to
                % flip them so that the vector generated from them is pointing
                % distally, like the vectors from the branches will be
                firstFive = fliplr(firstFive);
                % Convert the first and last points to coordinates
                [X,Y,Z] = ind2sub(size(aortaMask),firstFive);
                xu = X(end)-X(1);
                yu = Y(end)-Y(1);
                zu = Z(end)-Z(1);
                aortaUnitVector = 1/sqrt(xu^2 + yu^2 + zu^2) * [xu yu zu];
            end
            clear firstFive; clear X; clear Y; clear Z; clear xu; clear yu;
            clear zu;
        end
        clear originalNextNode; clear testNextNode;
        
        % Perform the scalar vector projection. In this case it is just the dot
        % product since we are projecting unit vectors.
        if(dot(aortaUnitVector,testNodeUnitVector) >= dot(aortaUnitVector,longestPathUnitVector))
            % Find the longestPathVector up to and including the commonNode
            idx = find(longestPathVector == commonNode);
            
            % Find the shortestPath from the common node (but not including it)
            % to the new end node
            idx2 = find(P == commonNode);
            
            % Rewrite the longest path
            longestPathVector = [longestPathVector(1:idx) fliplr(P(1:idx2-1))];
        end
    end
    %clear longestPathUnitVector; clear testNodeUnitVector; clear
    %aortaUnitVector; clear idx;
end

% Remove loops or wrinkles. Example: if the aorta goes from nodes:
%   A --> B --> C --> D --> E
%   But a link also exists from B --> D, change the longest path to:
%   A --> B --> D --> E
% Might be wise to calculate the true distance, and minimize distance
% instead of minimizing the number of nodes
longestPathVector = shortestpath(graph(A),longestPathVector(1),longestPathVector(end));

% Preallocate a structure to store the nodes and links that are important
%nodefinal = struct('idx',[],'links',[],'conn',[],'comx',[],'comy',[],'comz',[],'ep',[],'branch',[]);
%linkfinal = struct('n1',[],'n2',[],'point',[],'branch',[])

% Get the skeleton of the aorta
% Preallocate the array to store the skeleton
aortaskel = false(size(skel));
% Preallocate an array to store the indices (need to be in order!!)
aortaidx = zeros(sum(sum(sum(skel))),1);
nextindex = 1;
% Preallocate an array to store the aorta's links
aortalinks = zeros(length(longestPathVector)-1,1);
for ii = 1:(length(longestPathVector)-1)
    % Find the idx of the link to the next node
    node1 = longestPathVector(ii);
    node2 = longestPathVector(ii+1);
    idx = find(node(node1).conn == node2,1);
    aortaskel(link(node(node1).links(idx)).point) = true;
    if(link(node(node1).links(idx)).n1 == node1)
        % If the link is in the correct order
        aortaidx(nextindex:nextindex+length(link(node(node1).links(idx)).point)-1) = link(node(node1).links(idx)).point;
        nextindex = nextindex + length(link(node(node1).links(idx)).point);
        aortalinks(ii) = node(node1).links(idx);
    else
        % If the link is backwards
        aortaidx(nextindex:nextindex+length(link(node(node1).links(idx)).point)-1) = fliplr(link(node(node1).links(idx)).point);
        nextindex = nextindex + length(link(node(node1).links(idx)).point);
        % Store the link number as a negative if it needs to be reversed
        aortalinks(ii) = -1*node(node1).links(idx);
    end
    
end
clear node1; clear node2; clear idx; clear nextindex;
% Trim aortaidx (all indices are greater than zero)
aortaidx = aortaidx(aortaidx > 0);

%% Find the coordinates of the aorta's midline
% Preallocate an array to store the coordinates: X, Y, Z
orderedMidlineCoordinates = zeros(length(aortaidx),3);
[orderedMidlineCoordinates(:,1), orderedMidlineCoordinates(:,2), orderedMidlineCoordinates(:,3)] = ind2sub(size(aortaskel),aortaidx);
% Multiply by the voxelSize to get the coordinates instead of indicies
orderedMidlineCoordinates(:,1) = orderedMidlineCoordinates(:,1)*voxelSize(1);
orderedMidlineCoordinates(:,2) = orderedMidlineCoordinates(:,2)*voxelSize(2);
orderedMidlineCoordinates(:,3) = orderedMidlineCoordinates(:,3)*voxelSize(3);

% Make a subset of 12 points to trace a spline through (using too many
% points results in a wildly oscillating line due to the pixelation)
midlineSubset = orderedMidlineCoordinates(1:floor(size(orderedMidlineCoordinates,1)/11):end,:);
CS = cat(1,0,cumsum(sqrt(sum(diff(midlineSubset,[],1).^2,2))));
interpolatedMidline = interp1(CS, midlineSubset, unique([CS(:)' linspace(0,CS(end),200)]),'spline');

%% Find the three main branches of the aorta
% Find the branch points near the top of the aorta that are near the
% midline. These will have x values > the first node in longestPathVector,
% and will be a member of longestPathVector

tempidx = false(size(longestPathVector));
for ii = 2:length(tempidx)
    if(node(longestPathVector(ii)).comx < node(longestPathVector(1)).comx)
        % If the node is superior to the aortic root and on the aortic
        % midline (remember all the nodes in longestPathVector are on the
        % aortic midline
        tempidx(ii) = 1;
    end
end
possibleBranchNodes = longestPathVector(tempidx);
clear tempidx;

% Find the superior-most of the possible branch nodes
xmin = inf;
for ii = 1:length(possibleBranchNodes)
    if(node(possibleBranchNodes(ii)).comx < xmin)
        xmin = node(possibleBranchNodes(ii)).comx;
    end
end

% Get a list of endpoints that are superior to the highest of the branch
% nodes in possibleBranchNodes by at least one cm
tempidx = false(size(node));
% How many voxels is one cm?
onecm = 0;% ceil(10/voxelSize(1));
for ii = 1:length(tempidx)
    if(node(ii).comx < (xmin-onecm) && node(ii).ep == 1)
        % If the node is superior to the highest of the aorta branch nodes,
        % and is an endpoint
        tempidx(ii) = 1;
    end
end
possibleBranchEnds = find(tempidx);
clear tempidx; clear xmin; clear onecm;

% Look for and save the shortest path from each possible branch end to the
% aorta branch points
% Preallocate a structure to store the paths and distances
branchData = struct('path',{},'distance',{});
for ii = 1:length(possibleBranchEnds)
    % Allocate a structure to hold the temporary data
    searchData = struct('path',{},'distance',{});
    for jj = 1:length(possibleBranchNodes)
        [searchData(end+1).path,searchData(end+1).distance] = shortestpath(graph(A),possibleBranchEnds(ii),possibleBranchNodes(jj));
    end
    
    % Store the shortest path and distance for each potential endpoint
    [~,idx] = min(cat(1,searchData.distance));
    if(isempty(idx))
        branchData(ii).path = 0;
        branchData(ii).distance = inf;
    else
        branchData(ii).path = searchData(idx).path;
        branchData(ii).distance = searchData(idx).distance;
    end
end

% See which aorta branches are the sources of the three major branches
tempidx = false(max(possibleBranchNodes),1);
for ii = 1:size(branchData,2)
    if ~isempty(branchData(ii).path)
        tempidx(branchData(ii).path(end)) = true;
    end
end
branchNodes = find(tempidx > 0);
clear tempidx;

% Preallocate a vector to store the links associated with different
% branches (25 is a large enough number that it shouldn't be surpassed)
branchlinks = zeros(3,25);

% Make sure there are three branches
if(length(branchNodes) == 3)
    % Find which branch objects are part of which branch
    % First find which node is which branch, the lowest comy =
    % brachiocephalic, highest = left subclavian, middle = left common
    % carotid. Label 1 = brachiocephalic, 2 = left common carotid, 3 = left
    % subclavian. Is the code failing here? If so, consider comparing the
    % position of the node in the longestPathVector instead of the
    % coordinates
    temp = branchNodes;
    tempy = zeros(1,3);
    for ii = 1:3
        tempy(ii) = node(branchNodes(ii)).comy;
    end
    for ii = 1:3
        if(tempy(ii) == min(tempy))
            % Brachiocephalic
            temp(1) = branchNodes(ii);
        elseif(tempy(ii) == max(tempy))
            % Left subclavian
            temp(3) = branchNodes(ii);
        else
            % Left common carotid
            temp(2) = branchNodes(ii);
        end
    end
    branchNodes = temp;
    clear temp; clear tempy;
    
    % Copy aortaskel to a new variable
    labeledskel = int16(aortaskel);
    % For any link the the branches, give a label of 2
    for ii = 1:size(branchData,2)
        for jj = 1:length(branchData(ii).path)-1
            % Find the idx of the link to the next node
            node1 = min(branchData(ii).path(jj),branchData(ii).path(jj+1));
            node2 = max(branchData(ii).path(jj),branchData(ii).path(jj+1));
            idx = find(node(node1).conn == node2,1);
            % Label the links in the labeledSkel as 2
            if(branchData(ii).path(end) == branchNodes(1))
                % Brachiocephalic
                labeledskel(link(node(node1).links(idx)).point) = 2;
                % Write this link in the first zero element in row 1 of
                % branchlinks
                branchlinks(1,find(branchlinks(1,:) == 0,1)) = node(node1).links(idx);
            elseif(branchData(ii).path(end) == branchNodes(2))
                % Left common carotid
                labeledskel(link(node(node1).links(idx)).point) = 3;
                % Write this link in the first zero element in row 2 of
                % branchlinks
                branchlinks(2,find(branchlinks(2,:) == 0,1)) = node(node1).links(idx);
            else
                % Left subclavian
                labeledskel(link(node(node1).links(idx)).point) = 4;
                % Write this link in the first zero element in row 3 of
                % branchlinks
                branchlinks(3,find(branchlinks(3,:) == 0,1)) = node(node1).links(idx);
            end
        end
        
    end
    % elseif(length(branchNodes) == 2)
    % error('Only two branch nodes found, still need to implement')
    % Basically, should take the lowest comy as brachiocephalic
else
    % Attempt to extract just the brachiocephalic
    % Out of the potential branch ends, choose the one with the lowest
    % value of comy
    tempy = zeros(1,length(possibleBranchEnds));
    for ii = 1:length(possibleBranchEnds)
        tempy(ii) = node(possibleBranchEnds(ii)).comy;
    end
    % Find the branch end with the lowest y value
    [~,idx] = min(tempy);
    % Search for paths from this branch end to the aorta
    searchData = struct('path',{},'distance',{});
    for jj = 1:length(possibleBranchNodes)
        [searchData(end+1).path,searchData(end+1).distance] = shortestpath(graph(A),possibleBranchEnds(idx),possibleBranchNodes(jj));
    end
    clear idx;
    
    % Find the shortest path, this is the brachiocephalic (we assume)
    [~,baIdx] = min(cat(1,searchData.distance));
    branchData = struct('path',searchData(baIdx).path,'distance',searchData(baIdx).distance);
    
    % Copy aortaskel to a new variable
    labeledskel = int16(aortaskel);
    % For any link the the branches, give a label of 2
    for jj = 1:length(branchData(1).path)-1
        % Find the idx of the link to the next node
        node1 = min(branchData(1).path(jj),branchData(1).path(jj+1));
        node2 = max(branchData(1).path(jj),branchData(1).path(jj+1));
        idx = find(node(node1).conn == node2);
        % Label the links in the labeledSkel as 2
        % Brachiocephalic
        labeledskel(link(node(node1).links(idx)).point) = 2;
        % Write the link in the first row of branchlinks in the first
        % zero element
        branchlinks(1,find(branchlinks(1,:) == 0,1)) = node(node1).links(idx);
    end
end
clear tempy; clear baIdx;
% erase any columns of zeros
branchlinks(:,~any(branchlinks,1)) = [];

%% Generate the output datastruct
%struct.links...
%struct.nodes...
nodemap = NaN(size(node,2),2);
nodemap(:,1) = 1:size(nodemap,1);
linkmap = NaN(size(link,2),2);
linkmap(:,1) = 1:size(linkmap,1);
nodeindex = 1;
linkindex = 1;

% Find the new node number for the nodes in the aorta (IN ORDER)
for ii = 1:length(longestPathVector)
    nodemap(longestPathVector(ii),2) = nodeindex;
    nodeindex = nodeindex + 1;
end

% Find the new node number for the links in the aorta (IN ORDER)
for ii = 1:length(aortalinks)
    if(aortalinks(ii) > 0)
        linkmap(aortalinks(ii),2) = linkindex;
    else
        linkmap(abs(aortalinks(ii)),2) = -1*linkindex;
    end
    linkindex = linkindex + 1;
end

% Loop through the rest of the branch links
for ii = 1:3
    % For each aorta branch
    for jj = 1:size(branchlinks,2)
        if(branchlinks(ii,jj) > 0)
            linkmap(branchlinks(ii,jj),2) = linkindex;
            linkindex = linkindex + 1;
            % Add new node numbers if necessary
            % First node
            if(isnan(nodemap(link(branchlinks(ii,jj)).n1,2)))
                nodemap(link(branchlinks(ii,jj)).n1,2) = nodeindex;
                nodeindex = nodeindex + 1;
            end
            % Second node
            if(isnan(nodemap(link(branchlinks(ii,jj)).n2,2)))
                nodemap(link(branchlinks(ii,jj)).n2,2) = nodeindex;
                nodeindex = nodeindex + 1;
            end
        end
    end
end

% Clean up the maps/lookup tables
% Sort based on the new indicies
nodemapr = sortrows(nodemap,2);
% Eliminate rows with NaN (unused nodes)
nodemapr = nodemapr(~any(isnan(nodemapr),2),:);

% Eliminate rows with NaN (unused links)
linkmapr = linkmap(~any(isnan(linkmap),2),:);
% Sort based on the new indicies
% Add another column to store if reverse (ie if the number is negative).
% The next few lines ensure that the sorting is carried out normally while
% ignoring the negative signs if present (negative indicates the link is in
% the wrong direction and needs to be reversed)
temp = [linkmapr zeros(size(linkmapr,1),1)];
temp(:,3) = (temp(:,2) > 0);
temp(:,3) = temp(:,3) - (temp(:,2) < 0);
temp(:,2) = abs(temp(:,2));
temp = sortrows(temp,2);
temp(:,2) = temp(:,2) .* temp(:,3);
linkmapr = temp(:,1:2);
clear temp;

% Find which nodes are in which branches
% bcnodes stores the nodes in the brachiocephalic artery
% lccnodes stores the nodes in the left common carotid artery
% lsnodes stores the nodes in the left subclavian artery
bcnodes = [];
lccnodes = [];
lsnodes = [];
for ii = 1:size(branchData,2)
    if(branchData(ii).path(end) == branchNodes(1))
        % Brachiocephalic
        bcnodes = [bcnodes branchData(ii).path];
    elseif(branchData(ii).path(end) == branchNodes(2))
        % Left common carotid
        lccnodes = [lccnodes branchData(ii).path];
    elseif(branchData(ii).path(end) == branchNodes(3))
        % Left subclavian
        lsnodes = [lsnodes branchData(ii).path];
    else
        % Should not be reached
        warning('Unusual branch detected.')
    end
end
% Remove any duplicates
bcnodes = unique(bcnodes);
lccnodes = unique(lccnodes);
lsnodes = unique(lsnodes);

% Generate the new node structure
nodefinal = struct('idx',{},'links',{},'conn',{},'comx',{},'comy',{},'comz',{},'ep',{},'vessel',{});
for ii = 1:size(nodemapr,1)
    currentnode = node(nodemapr(ii,1));
    % idx unchanged
    % links
    templink = currentnode.links;
    for jj = 1:length(currentnode.links)
        % Rewrite the links to the new links
        templink(jj) = linkmap(currentnode.links(jj),2);
    end
    % eliminate NaN links
    currentnode.links = abs(templink(~isnan(templink) & templink ~= 0));
    % conn
    tempnode = currentnode.conn;
    for jj = 1:length(currentnode.conn)
        % Rewrite the nodes to the new nodes
        tempnode(jj) = nodemap(currentnode.conn(jj),2);
    end
    % eliminate NaN nodes
    currentnode.conn = tempnode(~isnan(tempnode) & tempnode ~= 0);
    % comx, comy, comz, ep unchanged
    % vessel
    oldnodenumber = nodemapr(ii,1);
    tempvessel = [];
    if(any(longestPathVector == oldnodenumber))
        tempvessel = 1;
    end
    if(any(bcnodes == oldnodenumber))
        tempvessel = [tempvessel 2];
    end
    if(any(lccnodes == oldnodenumber))
        tempvessel = [tempvessel 3];
    end
    if(any(lsnodes == oldnodenumber))
        tempvessel = [tempvessel 4];
    end
    currentnode.vessel = tempvessel;
    
    % Save the current node in the output
    nodefinal(end+1) = currentnode;
end

% Generate the new link structure
linkfinal = struct('n1',{},'n2',{},'point',{},'vessel',{});
for ii = 1:size(linkmapr,1)
    currentlink = link(linkmapr(ii,1));
    % n1: remap
    currentlink.n1 = nodemap(currentlink.n1,2);
    % n2: remap
    currentlink.n2 = nodemap(currentlink.n2,2);
    % point: same if linkmapr is positive, fliplr if negative
    if(linkmapr(ii,2) < 0)
        currentlink.point = fliplr(currentlink.point);
        % Switch n1 and n2
        temp = currentlink.n1;
        currentlink.n1 = currentlink.n2;
        currentlink.n2 = temp;
        clear temp;
    else
        % keep point the same, ie do nothing
    end
    % vessel
    oldlinknumber = linkmapr(ii,1);
    tempvessel = [];
    if(any(abs(aortalinks) == oldlinknumber))
        tempvessel = 1;
    end
    if(any(branchlinks(1,:) == oldlinknumber))
        tempvessel = [tempvessel 2];
    end
    if(any(branchlinks(2,:) == oldlinknumber))
        tempvessel = [tempvessel 3];
    end
    if(any(branchlinks(3,:) == oldlinknumber))
        tempvessel = [tempvessel 4];
    end
    currentlink.vessel = tempvessel;
    
    % Save the current link in the output
    linkfinal(end+1) = currentlink;
end

%% Plot the output if requested
if(plotflag)
    % Rerun the skeleton function to make [A, node, link] on
    % labeledskel(labeledskel > 0)
    [~,node2,link2] = Skel2Graph3D(labeledskel > 0,1);
    skel = Graph2Skel3D(node2,link2,w,l,h);
    
    % Generate a meshgrid
    [X,Y,Z] = meshgrid(voxelSize(2)*(1:size(aortaskel,2)),voxelSize(1)*(1:size(aortaskel,1)),voxelSize(3)*(1:size(aortaskel,3)));
    y = zeros(length(longestPathVector),1);
    x = zeros(length(longestPathVector),1);
    z = zeros(length(longestPathVector),1);
    for ii = 1:length(longestPathVector)
        y(ii) = node(longestPathVector(ii)).comy*voxelSize(2);
        x(ii) = node(longestPathVector(ii)).comx*voxelSize(1);
        z(ii) = node(longestPathVector(ii)).comz*voxelSize(3);
    end
    
    % Plot the skeleton
    figure();
    col=[.7 .7 .8];
    hiso = patch(isosurface(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
    hiso2 = patch(isocaps(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
    axis equal;axis off;
    lighting phong;
    isonormals(X,Y,Z,aortaMask,hiso);
    alpha(0.5);
    set(gca,'DataAspectRatio',[1 1 1])
    camlight;
    hold on;
    w=size(skel,1);
    l=size(skel,2);
    h=size(skel,3);
    [x,y,z]=ind2sub([w,l,h],find(labeledskel>0));
    plot3(voxelSize(2)*y,voxelSize(1)*x,voxelSize(3)*z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');
    set(gcf,'Color','white');
    view(140,80)
    
    [X,Y,Z] = meshgrid(voxelSize(2)*(1:size(aortaskel,2)),voxelSize(1)*(1:size(aortaskel,1)),voxelSize(3)*(1:size(aortaskel,3)));
    y = zeros(length(longestPathVector),1);
    x = zeros(length(longestPathVector),1);
    z = zeros(length(longestPathVector),1);
    for ii = 1:length(longestPathVector)
        y(ii) = node(longestPathVector(ii)).comy*voxelSize(2);
        x(ii) = node(longestPathVector(ii)).comx*voxelSize(1);
        z(ii) = node(longestPathVector(ii)).comz*voxelSize(3);
    end
    
    % Plot interpolated midline
    figure();
    col=[.7 .7 .8];
    hiso = patch(isosurface(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
    hiso2 = patch(isocaps(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
    axis equal;axis off;
    lighting phong;
    isonormals(X,Y,Z,aortaMask,hiso);
    alpha(0.5);
    set(gca,'DataAspectRatio',[1 1 1])
    camlight;
    hold on;
    plot3(y,x,z,'square','Markersize',8,'MarkerFaceColor','b','Color','k');
    plot3(interpolatedMidline(:,2),interpolatedMidline(:,1),interpolatedMidline(:,3),'LineWidth',2)
    set(gcf,'Color','white');
    view(140,80)
    
    
    
    % Plot the skeleton with branches labeled
    figure();
    col=[.7 .7 .8];
    hiso = patch(isosurface(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
    hiso2 = patch(isocaps(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
    axis equal;axis off;
    lighting phong;
    isonormals(X,Y,Z,aortaMask,hiso);
    alpha(0.5);
    set(gca,'DataAspectRatio',[1 1 1])
    camlight;
    hold on;
    w=size(skel,1);
    l=size(skel,2);
    h=size(skel,3);
    if(max(max(max(labeledskel)))==4)
        [xa,ya,za]=ind2sub([w,l,h],find(labeledskel == 1));
        [xb,yb,zb]=ind2sub([w,l,h],find(labeledskel == 2));
        [xcc,ycc,zcc]=ind2sub([w,l,h],find(labeledskel == 3));
        [xls,yls,zls]=ind2sub([w,l,h],find(labeledskel == 4));
        p1 = plot3(voxelSize(2)*ya,voxelSize(1)*xa,voxelSize(3)*za,'square','Markersize',4,'MarkerFaceColor','r','Color','r');
        p2 = plot3(voxelSize(2)*yb,voxelSize(1)*xb,voxelSize(3)*zb,'square','Markersize',4,'MarkerFaceColor','k','Color','k');
        p3 = plot3(voxelSize(2)*ycc,voxelSize(1)*xcc,voxelSize(3)*zcc,'square','Markersize',4,'MarkerFaceColor','b','Color','b');
        p4 = plot3(voxelSize(2)*yls,voxelSize(1)*xls,voxelSize(3)*zls,'square','Markersize',4,'MarkerFaceColor','g','Color','g');
        legend([p1 p2 p3 p4],'Aorta','Brachiocephalic','Left Common Carotid','Left subclavian')
    else
        [xa,ya,za]=ind2sub([w,l,h],find(labeledskel == 1));
        [xb,yb,zb]=ind2sub([w,l,h],find(labeledskel == 2));
        p1 = plot3(voxelSize(2)*ya,voxelSize(1)*xa,voxelSize(3)*za,'square','Markersize',4,'MarkerFaceColor','r','Color','r');
        p2 = plot3(voxelSize(2)*yb,voxelSize(1)*xb,voxelSize(3)*zb,'square','Markersize',4,'MarkerFaceColor','k','Color','k');
        legend([p1 p2],'Aorta','Brachiocephalic')
    end
    set(gcf,'Color','white');
    view(140,80)
end
end