function cardiac4dflow(varargin)
% main.m
%
%   This code takes in a NUNDAin.mat struct and the scans on NUNDA and does
%   some automated analysis. It fits a midline and planes, and extracts the
%   ascending aorta. Then it uses the planes to find hemodynamic
%   parameters. It also uses a user-supplied (within NUNDAin.mat) mask to
%   find global hemodynamic and volumetric parameters. Results are saved
%   either in the selected folder, or in the user supplied path. 
%
%   Inputs: see example usage, below....
%   - inputpath (input1): patient folder. The code looks at path/NUNDA to 
%   find the file NUNDAin.mat. This structure should be generated using
%   generateNUNDAinput.m. If not provided, user will need to use the GUI to
%   select the input path.
%   - outputpath (input2): folder for storing results. If not provided, the
%   results will be stored in a timestamped folder within the inputpath.
%   - list of 4d flow scan ids
%
%   Constants:
%   - numberOfPlanes: the number of planes placed in the plane analysis
%   - spacing: the isotropic spacing on the plane, in mm
%
%   Outputs:
%   - NUNDAout.mat: a structure containing the results of the calculations
%
%   Example usage:
%   - main(): with no inputs the user needs to select the patient folder
%   - main(path): the code will look within the path to find a 'NUNDA' 
%   folder containing NUNDAin.mat. Results will be saved in a subfolder
%   created in the path folder
%   - main(inputpath,outputpath): the code will look at inputpath/NUNDA for
%   the NUNDAin.mat file, and save outputs in the output folder supplied.
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
%   - Need to add a field in the NUNDAout structure that gives the units
%   for the fields in the output
%   - Need to have to code also use CEMRA data if provided
%   - Check that the hemodynamic parameters (particularly retrograde flow,
%   cardiac output, and stroke volume are reasonable)
%   - Solve skeletonization issues with problem patients.
%
% Written by Mike Scott, August 2017 (Northwestern University)
% michael.scott1@northwestern.edu

fprintf('==================================================\n');
fprintf('===              Aneurysm Analysis             ===\n');
fprintf('==================================================\n');

addpath(genpath([pwd() filesep() 'functions']));
formatOut = 'yyyymmdd_HHMMSS';
if nargin == 0
    folder_name = uigetdir('C:\Users\Mike\Desktop\Working','Select the patient folder:');
    results_folder = [folder_name filesep() datestr(now,formatOut) '_results'];
elseif nargin == 1
    folder_name = varargin{1};
    results_folder = [folder_name filesep() datestr(now,formatOut) '_results'];
elseif nargin == 2
    folder_name = varargin{1};
    results_folder = varargin{2};
elseif nargin == 3
    folder_name = varargin{1};
    results_folder = varargin{2};
    folderlistcell = varargin{3};
else
    error('Too many input arguments!')
end

% Get the input
NUNDApath = [folder_name filesep() 'NUNDA' filesep() 'NUNDAin.mat'];

% Create a folder to store the results
if ~exist(results_folder,'dir')
    mkdir(results_folder)
end

fprintf(' - Loading data from:\n   %s\n',folder_name);
tic;

if(exist(NUNDApath,'file') == 2)
    % Load the NUNDA input
    load(NUNDApath);
    if(isfield(NUNDAin,'ceMask'))
        % Implement later
    end
    if(isfield(NUNDAin,'pcMask'))
        velStruct = NUNDAin.velStruct;
        magStruct = NUNDAin.magStruct;
        mask_struct = NUNDAin.pcMask;
        aortaMask = NUNDAin.pcMask.mask;
        vox = mask_struct.vox;
        % Need for now
        isCEMRA = false;
    end
else
    error('No NUNDA input exists in the chosen folder (check for a NUNDA subfolder...)')
end

% Try to find the date and nominal interval
if exist('folderlistcell','var')
    idx4Dflow = 1:length(folderlistcell);
else
    folderlist = dir(folder_name);
    folderlistcell = {folderlist(:).name};
    idx4D = contains(folderlistcell,'4D','IgnoreCase',true);
    idxflow = contains(folderlistcell,'flow','IgnoreCase',true);
    % flowand4D represents folders that have 4D and flow in their name
    flowand4D = idx4D | idxflow;
    idx4Dflow = find(flowand4D == 1);
end

% Search for the scan date and nominal interval
for ii = 1:length(idx4Dflow)
    % Get a list of dicom files in the folder
    dicomlist = dir([folder_name filesep() folderlistcell{idx4Dflow(ii)} filesep() '*.dcm']);
    info = dicominfo([folder_name filesep() folderlistcell{idx4Dflow(ii)} filesep() dicomlist(1).name]);
    
    if(isfield(info,'NominalInterval') && info.NominalInterval ~= 0)
        nominalInterval = info.NominalInterval;
    end
    if(isfield(info,'StudyDate'))
       studyDate = info.StudyDate; 
    end
    if(exist('studyDate','var') && exist('nominalInterval','var'))
       break; 
    end
end
clear folderlist; clear folderlistcell; clear idx4D; clear idxflow; 
clear idx4Dflow;


% Ensure that aortamask has only one volume
[L, NUM] = bwlabeln(aortaMask);
if(NUM > 1)
    % If there are two regions, assume one is the three branches and the
    % descending aorta (the larger), and the smaller one is the missing
    % piece of the ascending aorta.
    tempproperties = regionprops(L);
    [~,idx] = max(cat(1,tempproperties.Area));
    aortaMask = (L==idx);
end
clear tempproperties; clear L; clear NUM; clear idx;

% Apply the aortaMask to the velocity data if the data is velStruct
if ~isCEMRA
    velData = velStruct.dataAy .* repmat(aortaMask, [1 1 1 3 size(velStruct.dataAy,5)]); 
end

% Run the skeletonization scheme
try
    [labeledskel, interpolatedMidline, node, link] = skeletonize_Aorta_NUNDA(mask_struct);
    fprintf(' - Aorta skeletonization routine completed.\n');
catch
    % If the skeletonization fails, erode the mask struct and try again
    fprintf(' - Aorta skeletonization routine failed, retrying after erosion.\n');
    % Erode the mask
    se = strel(ones(3,3,3));
    erodedMask = imerode(mask_struct.mask,se);
    % Get rid of small islands that were created
    CC = bwconncomp(erodedMask);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    BW = false(size(erodedMask));
    BW(CC.PixelIdxList{idx}) = 1;
    erodedMask = BW;
    
    originalMask = mask_struct.mask;
    mask_struct.mask = erodedMask;
    % Retry the skeletonization
    [labeledskel, interpolatedMidline, node, link] = skeletonize_Aorta_NUNDA(mask_struct);
    % Restore the mask
    mask_struct.mask = originalMask;
    fprintf(' - Aorta skeletonization routine completed after erosion.\n');
end

% Find the index of the brachiocephalic junction
for ii = 1:size(node,2)
    if(length(node(ii).vessel) == 2 && node(ii).vessel(1) == 1 && node(ii).vessel(2) == 2)
        junctionIdx = ii;
    end
end

% Find the coordinates of the junction (P0), as well as the points along
% the aorta immediately preceding (P1) and after the junction (P2). Instead
% of using the exact point, find the closest point on interpolated midline!
junctionCoords = [node(junctionIdx).comx*vox(1) node(junctionIdx).comy*vox(2) node(junctionIdx).comz*vox(3)];
% Calculate the distance from the junction to all the points in the midline
distances = sqrt(sum(bsxfun(@minus, interpolatedMidline, junctionCoords).^2,2));
P0 = interpolatedMidline(distances==min(distances),:);
[~,P0idx] = min(distances);


% Find the plane that is perpendicular to the midline at the junction point
P1 = interpolatedMidline(P0idx-1,:);
P2 = interpolatedMidline(P0idx+1,:);
P10 = P1-P0;
P20 = P2-P0;
N = dot(P10,P10)*P20-dot(P20,P20)*P10; % <-- Approx. tangent direction
Nhat = N ./ sqrt(sum(N.^2));
%T = null(N).'; % Get two orthogonal unit vectors which are orthog. to N

% Find the equation of the plane: Nhat is the unit normal, expressed as
%   nhat = (a,b,c)
% The plane has equation a*x + b*y + c*z + d = 0
% Need to solve for d = -a*x0 - b*y0 - c*z0, where (x0,y0,z0) is a point on
% the plane (use the point P0).
d = -Nhat(1) * P0(1) - Nhat(2) * P0(2) - Nhat(3) * P0(3);
planeVec = [Nhat d];
pointVec = [0 0 0 1];

% Method:
% 1) Check if the sinus of valsalva has a positive value when the
% coordinates are dot producted with planevec (note the coordinates need to
% be expressed as (x,y,z,1)
% 2) Check all points (or all aorta points?), all points on the same side
% of the plane as the sinus of valsalva will have the same sign

% Sinus of Valsalva is the first point on interpolatedMidline
pointVec(1:3) = interpolatedMidline(1,:);
signCheck = dot(pointVec,planeVec);

% Make a meshgrid of the mm coordinates of the image
x = (1:size(labeledskel,1))*vox(1);
y = (1:size(labeledskel,2))*vox(2);
z = (1:size(labeledskel,3))*vox(3);
[Y,X,Z] = meshgrid(y,x,z);
clear x; clear y; clear z;

% Preallocate an array to store the mask
planeMask = false(size(labeledskel));
maxIndex = size(labeledskel,1)*size(labeledskel,2)*size(labeledskel,3);
if(signCheck > 0)
    % If the sign of the dot product of the normal vector to the plane and
    % the point at the end of the midline in the SOV is positive, all
    % points with positive dot products with the normal vector will be on
    % the same side of the plane as the SOV
    for ii = 1:maxIndex
        pointVec(1:3) = [X(ii) Y(ii) Z(ii)];
        if(dot(pointVec,planeVec) > 0)
            planeMask(ii) = 1;
        end
    end
elseif(signCheck < 0)
    % If the sign of the dot product of the normal vector to the plane and
    % the point at the end of the midline in the SOV is negative, all
    % points with negative dot products with the normal vector will be on
    % the same side of the plane as the SOV
    for ii = 1:maxIndex
        pointVec(1:3) = [X(ii) Y(ii) Z(ii)];
        if(dot(pointVec,planeVec) < 0)
            planeMask(ii) = 1;
        end
    end
else
    error('Sinus of Valsalva is on the same plane as the brachiocephalic junction')
end

% Apply the mask
maskedVolume = planeMask & aortaMask;
% Use bwlabeln to find the number of objects
[L, ~] = bwlabeln(maskedVolume);
% Find the labeled region which includes the sinus of Valsalva (Node 1 of
% the structure node)
regionToKeep = L(node(1).comx,node(1).comy,node(1).comz);
maskedVolume = (L == regionToKeep);

% Add back in any area that was deleted when the plane cut the ascending
% aorta (ie if the plane cut a portion of the sinus off)
% Subtract the maskedVolume (AA only) from the aortaMask (full mask)
tempVol = aortaMask - maskedVolume;
% Label the volumes
[L, NUM] = bwlabeln(tempVol);
if(NUM == 2)
    % If there are two regions, assume one is the three branches and the
    % descending aorta (the larger), and the smaller one is the missing
    % piece of the ascending aorta.
    tempproperties = regionprops(L);
    [~,idx] = min(cat(1,tempproperties.Area));
    tempVol = (L==idx);
    maskedVolume = maskedVolume | tempVol;
elseif(NUM>2)
    tempproperties = regionprops(L);
    maxVol = max(cat(1,tempproperties.Area));
    for ii = 1:NUM
        if tempproperties(ii).Area ~= maxVol
           % Ignore the largest piece, almost certainly the descending 
           % aorta, then check each smaller piece. Add it back to the main
           % chunk of the ascending aorta and check connectivity, if only
           % one region is found, keep the addition, otherwise, ignore it
           % Add the region to the masked volume (bulk of ascending aorta)
           tempVol = maskedVolume | (L==ii);
           % Count the number of regions
           [~, NUM_check] = bwlabeln(tempVol);
           if NUM_check == 1
               % If the two regions added together form one region, keep it
               maskedVolume = tempVol;
           end
        end
    end
    
end
clear tempVol; clear tempproperties; clear idx; clear maxVol;
fprintf(' - Mask applied to isolate the ascending aorta.\n');



% Get properties of the ascending aorta
AAproperties = regionprops(maskedVolume);
% Find the volume in mm^3
AAvolume = vox(1) * vox(2) * vox(3) * AAproperties.Area;
fprintf(' - Ascending aorta volume calculated: %g mL.\n',AAvolume/1000);


%% Do the diameter calculation
if exist('mrstruct_mask', 'var')
    maskedVolume_mask = mrstruct_mask;
    maskedVolume_mask.dataAy = maskedVolume;
else
    maskedVolume_mask = mask_struct;
    maskedVolume_mask.dataAy = maskedVolume;
end

% Use evalc to suppress output from pimDiameterMap
fprintf(' - Starting parallel pool for diameter map calculation...\n');
[~,Diameter_matrix] = evalc('pimDiameterMap(maskedVolume_mask,1,1,0,0)');
fprintf(' - Diameter map calculation complete.\n');

% Generate a mask that will eliminate the end caps, since they distort the
% volume measurements
% The erosion of the half-plane mask (planeMask) will take off the pixels
% at the distal end of the ascending aorta.
se = strel(ones(3,3,3));
erodedMask = imerode(planeMask,se);
% Generate a half plane mask that cuts off the proximal end of the aortic 
% root
% Find the plane that is perpendicular to the midline at the junction point
P0 = interpolatedMidline(2,:);
P1 = interpolatedMidline(1,:);
P2 = interpolatedMidline(3,:);
P10 = P1-P0;
P20 = P2-P0;
N = dot(P10,P10)*P20-dot(P20,P20)*P10; % <-- Approx. tangent direction
Nhat = N ./ sqrt(sum(N.^2));
T = null(N).'; % Get two orthogonal unit vectors which are orthog. to N

% Find the equation of the plane: Nhat is the unit normal, expressed as
%   nhat = (a,b,c)
% The plane has equation a*x + b*y + c*z + d = 0
% Need to solve for d = -a*x0 - b*y0 - c*z0, where (x0,y0,z0) is a point on
% the plane (use the point P0).
d = -Nhat(1) * P0(1) - Nhat(2) * P0(2) - Nhat(3) * P0(3);
planeVec = [Nhat d];
pointVec = [0 0 0 1];

% The third point on the midline is in the half-space that should be kept
pointVec(1:3) = interpolatedMidline(3,:);
signCheck = dot(pointVec,planeVec);

% Preallocate an array to store the mask
planeMask2 = false(size(labeledskel));
maxIndex = size(labeledskel,1)*size(labeledskel,2)*size(labeledskel,3);
if(signCheck > 0)
    % If the sign of the dot product of the normal vector to the plane and
    % the point at the end of the midline in the SOV is positive, all
    % points with positive dot products with the normal vector will be on
    % the same side of the plane as the SOV
    for ii = 1:maxIndex
        pointVec(1:3) = [X(ii) Y(ii) Z(ii)];
        if(dot(pointVec,planeVec) > 0)
            planeMask2(ii) = 1;
        end
    end
elseif(signCheck < 0)
    % If the sign of the dot product of the normal vector to the plane and
    % the point at the end of the midline in the SOV is negative, all
    % points with negative dot products with the normal vector will be on
    % the same side of the plane as the SOV
    for ii = 1:maxIndex
        pointVec(1:3) = [X(ii) Y(ii) Z(ii)];
        if(dot(pointVec,planeVec) < 0)
            planeMask2(ii) = 1;
        end
    end
else
    error('Check the code, point 3 of the midline should not be on the plane')
end

% Combine the masks
eraserMask = erodedMask & planeMask2;
% Apply the masks to Diameter_matrix
Diameter_matrix = Diameter_matrix .* (planeMask2 & erodedMask);

fprintf(' - Ends of the diameter mask clipped (they introduce error).\n');

% For each point in the Diameter_matrix, find the closest point in
% interpolatedMidline. Remember to work in mm coordinates, since the voxels
% can be anisotropic. Store the data in a structure
% D_struct
% .idx: idx in the Diameter_matrix
% .diameter (in cm): from Pim's code
% .midlineidx: idx of nearest midline point in interpolatedMidline
% .radius (in cm): estimated from calculating the distance to the midline
% Note that radius will be an overestimate, because it is calculating the
% distance to the nearest midline point, which is likely not perpendicular.
% Still, assuming the radius is large compared to the spacing between the
% midline points, this should be a reasonable estimate.

D_struct = struct('idx',{},'diameter',{},'midlineIdx',{},'radius',{});
% Get the idx of nonzero elements of Diameter_matrix
diameterIdx = find(Diameter_matrix > 0);
% Generate a subset of the interpolatedMidline data (to speed up the
% search)
midlineSubset = interpolatedMidline(1:P0idx,:);
% temp = midlineSubset(:,2);
% midlineSubset(:,2) = midlineSubset(:,1);
% midlineSubset(:,1) = temp;
distances = zeros(size(midlineSubset,1),1);
diameterCoords = [X(diameterIdx) Y(diameterIdx) Z(diameterIdx)];
for ii = 1:length(diameterIdx)
   % Store the idx
   D_struct(ii).idx = diameterIdx(ii);
   % Get the coordinates
   currentCoords = diameterCoords(ii,:);
   % Store the diameter calculated in Pim's code
   D_struct(ii).diameter = Diameter_matrix(diameterIdx(ii));
   % Calculate distances through P0idx (the last point kept in the
   % ascending aorta)
   for jj = 1:size(midlineSubset,1)
       distances(jj) = sqrt((currentCoords(1) - midlineSubset(jj,1))^2 + (currentCoords(2) - midlineSubset(jj,2))^2 + (currentCoords(3) - midlineSubset(jj,3))^2)/10; % /10 for cm
   end
   % Find the minimim distance and the index, store both
   [D_struct(ii).radius, D_struct(ii).midlineIdx] = min(distances);
end

% Bin the data by midline point
maxPoint = max(cat(1,D_struct.midlineIdx));
plotData = struct('idx',{},'diameters',{},'maxD',{},'minD',{},'averageD',{},'distance',{},'percentDistance',{});
diamData = cat(1,D_struct.diameter);
idxData = cat(1,D_struct.midlineIdx);
for ii = 1:maxPoint
   plotData(ii).idx = ii;
   plotData(ii).diameters = diamData(idxData == ii);
   plotData(ii).maxD = max(plotData(ii).diameters);
   if(isempty(plotData(ii).maxD))
       plotData(ii).maxD = NaN;
   end
   plotData(ii).minD = min(plotData(ii).diameters);
   if(isempty(plotData(ii).minD))
       plotData(ii).minD = NaN;
   end
   plotData(ii).averageD = mean(plotData(ii).diameters);
   if ii == 1
       plotData(ii).distance = 0;
   else
       addDistance = sqrt((midlineSubset(ii,1) - midlineSubset(ii-1,1))^2 + (midlineSubset(ii,2) - midlineSubset(ii-1,2))^2 + (midlineSubset(ii,3) - midlineSubset(ii-1,3))^2);
       plotData(ii).distance = plotData(ii-1).distance + addDistance;
   end
end
for ii = 1:size(plotData,2)
   plotData(ii).percentDistance = plotData(ii).distance / plotData(end).distance * 100;
end
fprintf(' - Diameter map data binned to nearest point on the midline.\n');

%% Plane analysis
numberOfPlanes = 10;
spacing = 1; % mm
fprintf(' - Beginning plane analysis with %i planes\n    and %d mm grid spacing for interpolation.\n',numberOfPlanes,spacing);
% Place the planes on points in the midline subset
planeMidlineIdx = round(linspace(2,size(midlineSubset,1)-1,numberOfPlanes));
planeStruct = struct('coords',[],'xvec',[],'yvec',[],'zvec',[],'Nhat',[],'Vx',[],'Vy',[],'Vz',[],'Vparallel',[],'MajorAxis',[],'MinorAxis',[]);
for ii = 1:numberOfPlanes
    % Get the maximum diameter at that midlineIdx
    % *10 for mm, *1.1 to ensure that the plane covers the entire area
    if planeMidlineIdx(ii) == 1
        maxWidth = ceil(plotData(2).maxD*10*1.10);
    elseif planeMidlineIdx(ii) <= size(plotData,2)
        maxWidth = ceil(plotData(planeMidlineIdx(ii)).maxD*10*1.10);
        if isnan(maxWidth)
           maxWidth = 0; 
        end
    else
        maxWidth = ceil(plotData(end).maxD*10*1.10);
    end
    % Get the normal vector at that location
    % Find the plane that is perpendicular to the midline at the junction point
    P0 = midlineSubset(planeMidlineIdx(ii),:);
    P1 = midlineSubset(planeMidlineIdx(ii)-1,:);
    P2 = midlineSubset(planeMidlineIdx(ii)+1,:);
    P10 = P1-P0;
    P20 = P2-P0;
    N = dot(P10,P10)*P20-dot(P20,P20)*P10; % <-- Approx. tangent direction
    Nhat = N ./ sqrt(sum(N.^2));
    
    % Input the coordinates, normal vector, plane size, and spacing
    planeCoordinates = findPlaneCoordinates(P0, Nhat, maxWidth, spacing);
    planeStruct(ii).coords = planeCoordinates;
    planeStruct(ii).Nhat = Nhat;
    
    % Convert the query points into x,y,z
    xvec = zeros(size(planeCoordinates,1)*size(planeCoordinates,2),1);
    yvec = xvec;
    zvec = xvec;
    
    index = 1;
    for jj = 1:size(planeCoordinates,2)
        for kk = 1:size(planeCoordinates,1)
            xvec(index) = planeCoordinates(kk,jj,1);
            yvec(index) = planeCoordinates(kk,jj,2);
            zvec(index) = planeCoordinates(kk,jj,3);
            index = index + 1;
        end
    end
    planeStruct(ii).xvec = xvec;
    planeStruct(ii).yvec = yvec;
    planeStruct(ii).zvec = zvec;
end

fprintf(' - Plane coordinates calculated, beginning interpolation...\n');

if ~isCEMRA
    % Generate a meshgrid
    [X,Y,Z] = meshgrid(vox(2)*(1:size(labeledskel,2)),vox(1)*(1:size(labeledskel,1)),vox(3)*(1:size(labeledskel,3)));
    %F = griddedInterpolant(X,Y,Z,V);
    for ii = 1:numberOfPlanes
        Nhat = planeStruct(ii).Nhat;
        % Preallocate arrays to store the interpolated data
        Vx = zeros(size(planeStruct(ii).coords,1),size(planeStruct(ii).coords,2),size(velStruct.dataAy,5));
        Vy = zeros(size(planeStruct(ii).coords,1),size(planeStruct(ii).coords,2),size(velStruct.dataAy,5));
        Vz = zeros(size(planeStruct(ii).coords,1),size(planeStruct(ii).coords,2),size(velStruct.dataAy,5));
        Vp = zeros(size(planeStruct(ii).coords,1),size(planeStruct(ii).coords,2),size(velStruct.dataAy,5));
        
        % Generate n to simplify code on the next few lines
        n = size(Vx,1);
           
        % Interpolate
        for tt = 1:size(velStruct.dataAy,5)
            Vx(:,:,tt) = interpolateVelocity(X,Y,Z,velStruct.dataAy(:,:,:,1,tt).*aortaMask,planeStruct(ii).coords);
            Vy(:,:,tt) = interpolateVelocity(X,Y,Z,velStruct.dataAy(:,:,:,2,tt).*aortaMask,planeStruct(ii).coords);
            Vz(:,:,tt) = interpolateVelocity(X,Y,Z,velStruct.dataAy(:,:,:,3,tt).*aortaMask,planeStruct(ii).coords);
            
            % Calculate vparallel (note Nhat updated in the outer loop!)
            Vp(:,:,tt) = dot(cat(3,Vx(:,:,tt),Vy(:,:,tt),Vz(:,:,tt)),cat(3,repmat(Nhat(1),n,n),repmat(Nhat(2),n,n),repmat(Nhat(3),n,n)),3);
        end
        
        % Ensure that the plane is only one region (take the largest)
        Vx = repmat(bwareafilt((Vx(:,:,1) ~= 0),1),1,1,size(Vx,3)).*Vx;
        Vy = repmat(bwareafilt((Vx(:,:,1) ~= 0),1),1,1,size(Vy,3)).*Vy;
        Vz = repmat(bwareafilt((Vx(:,:,1) ~= 0),1),1,1,size(Vz,3)).*Vz;
        Vp = repmat(bwareafilt((Vx(:,:,1) ~= 0),1),1,1,size(Vp,3)).*Vp;
        
        
        % Set the zero values to NaN to shape the planes
        Vx(Vx == 0) = NaN;
        Vy(Vy == 0) = NaN;
        Vz(Vz == 0) = NaN;
        Vp(Vp == 0) = NaN;
        
        % Store the data
        planeStruct(ii).Vx = Vx;
        planeStruct(ii).Vy = Vy;
        planeStruct(ii).Vz = Vz;
        
        planeStruct(ii).Vparallel = Vp;
        
        % Calculate the major and minor axis using one of the velocity
        % components
        planeProps = regionprops(~isnan(planeStruct(ii).Vx(:,:,1)),'MajorAxisLength','MinorAxisLength');
        if size(planeProps,1) > 1
            % Assume the object with the longest major axis is the correct
            % one
            [~,idx] = max(cat(1,planeProps.MajorAxisLength));
            planeStruct(ii).MajorAxis = spacing * planeProps(idx).MajorAxisLength;
            planeStruct(ii).MinorAxis = spacing * planeProps(idx).MinorAxisLength;
        else
            planeStruct(ii).MajorAxis = spacing * planeProps.MajorAxisLength;
            planeStruct(ii).MinorAxis = spacing * planeProps.MinorAxisLength;
        end
    end
end


%% Calculate the MIP
if ~isCEMRA
    [veloMIP,vmaxcoords,tsystole] = velocityMIP(velStruct.dataAy,aortaMask);
    % Compute the MIP for the ascending aorta
    [AAveloMIP,AAvmaxcoords] = velocityMIP(velStruct.dataAy,maskedVolume);
    
    % Plot the MIP
    % Tile the MIPs (note that MIPx is oriented 90 degrees rotated)
    MIPplot = NaN(size(veloMIP.MIPz,1) + size(veloMIP.MIPx,2),size(veloMIP.MIPz,2) + size(veloMIP.MIPy,2));
    MIPplot(1:size(veloMIP.MIPz,1),1:size(veloMIP.MIPz,2)) = veloMIP.MIPz;
    MIPplot(size(veloMIP.MIPz,1)+1:end,1:size(veloMIP.MIPz,2)) = veloMIP.MIPx';
    MIPplot(1:size(veloMIP.MIPz,1),size(veloMIP.MIPz,2)+1:end) = veloMIP.MIPy;
    
    % Find the coordinates for the MIP
    MRsize = size(velStruct.dataAy);
    
    % Make the plotting vectors for the MIP (note uneven spacing due to
    % voxel anisotropy)
    xmips = [linspace(vox(1),vox(1)*MRsize(1),MRsize(1)) (linspace(vox(3),vox(3)*MRsize(3),MRsize(3)) + vox(1)*MRsize(1))];
    ymips = [linspace(vox(2),vox(2)*MRsize(2),MRsize(2)) (linspace(vox(3),vox(3)*MRsize(3),MRsize(3)) + vox(2)*MRsize(2))];
    zmips = zeros(length(xmips),length(ymips));
    
    % Find the coordinates for the point of maximum velocity
    coords1 = veloMIP.MIPzcoords .* [vox(1) vox(2)];
    %coords2(1) = vox(1)*(MRsize(1)-1) + veloMIP.MIPxcoords(2) * vox(3);
    coords2(1) = vox(1)*(MRsize(1)) + veloMIP.MIPxcoords(2) * vox(3);
    coords2(2) = coords1(2);
    coords3(1) = coords1(1);
    %coords3(2) = vox(2)*(MRsize(2)-1) + veloMIP.MIPycoords(2) * vox(3);
    coords3(2) = vox(2)*(MRsize(2)) + veloMIP.MIPycoords(2) * vox(3);
    yplot = [coords1(1) coords2(1) coords3(1)];
    xplot = [coords1(2) coords2(2) coords3(2)];
    
    % Plot the MIP
    fig1 = figure(1);
    set(gcf,'color','k');
    title('title')
    surf(ymips,xmips,zmips,MIPplot,'EdgeColor','none','FaceColor','interp')
    set(gca,'Ydir','reverse')
    colormap('jet')
    c = colorbar;
    c.Label.String = 'Velocity (m/s)';
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'Visible','off');    
    set(c,'color','white')
    set(c,'Limits',[0 max(max(MIPplot))])
    view(2)
    hold on;
    plot(xplot,yplot,'ow')
    daspect([1 1 1])
    str = sprintf('V_{max} %.3f m/s',veloMIP.MIPz(veloMIP.MIPzcoords(1),veloMIP.MIPzcoords(2)));
    text(ymips(end),xmips(end),str,'Color','white','FontSize',14,'Interpreter','tex','HorizontalAlignment','right','VerticalAlignment','bottom');
    str = sprintf('Aorta Velocity MIP: %s',[velStruct.patient]);
    text(ymips(1),xmips(1),str,'Color','white','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
    % Save the figure
    saveas(fig1,[results_folder filesep() 'veloMIP.png'])
end

%% Plane analysis
if ~isCEMRA
    if(exist('nominalInterval','var'))
        planeStruct = planeAnalysis(planeStruct,nominalInterval);
    else
        planeStruct = planeAnalysis(planeStruct,velStruct.tr);
    end
end

%% Generate a structure containing data to be saved
NUNDAout.patient = NUNDAin.patient;
NUNDAout.planes = planeStruct;
NUNDAout.AAVmax = AAveloMIP.MIPz(AAveloMIP.MIPzcoords(1),AAveloMIP.MIPzcoords(2));
NUNDAout.AAVmaxidx = AAvmaxcoords;
NUNDAout.AAvolume = AAvolume;
NUNDAout.Vmax = veloMIP.MIPz(veloMIP.MIPzcoords(1),veloMIP.MIPzcoords(2));
NUNDAout.Vmaxidx = vmaxcoords;
NUNDAout.vox = vox;
NUNDAout.patientPath = folder_name;
if(exist('nominalInterval','var'))
    NUNDAout.NominalInterval = nominalInterval;
end
NUNDAout.StudyDate = studyDate;
% Save the structure
save([results_folder filesep() 'NUNDAout.mat'],'NUNDAout');

%% Temporary plotting
fprintf(' - Plotting results...\n\n');
%% Plot the diameter data (simple)
fig2 = figure(2);
%clf(fig2)
plot(cat(1,plotData.percentDistance),cat(1,plotData.maxD),'-r')
hold on
plot(cat(1,plotData.percentDistance),cat(1,plotData.averageD),'-g')
plot(cat(1,plotData.percentDistance),cat(1,plotData.minD),'-b')
legend('Max diameter', 'Mean diameter', 'Min diameter')
legend('Location','northwest')
axis([0 100 0 max(cat(1,plotData.maxD))])
xlabel('Percent distance along ascending aorta (%)')
ylabel('Diameter (cm)')

% Save the figure
saveas(fig2,[results_folder filesep() 'diam.png'])

%% Diameters Boxplot
xx = cat(1,D_struct.midlineIdx);
yy = cat(1,D_struct.diameter);
fig3 = figure(3);
boxplot(yy,xx)
title('Diameter distribution by midline point')
xlabel('Midline point')
ylabel('Diameter (cm)')
% Save the figure
saveas(fig3,[results_folder filesep() 'diamBoxplot.png'])

%% Plot the skeleton
% Rerun the skeleton function to make [A, node, link] on
% labeledskel(labeledskel > 0)
% Eliminate unconnected nodes
w = size(labeledskel,1);
l = size(labeledskel,2);
h = size(labeledskel,3);
[~,node2,link2] = Skel2Graph3D(labeledskel > 0,1);
skel = Graph2Skel3D(node2,link2,w,l,h);

% Generate a meshgrid
[X,Y,Z] = meshgrid(vox(2)*(1:size(skel,2)),vox(1)*(1:size(skel,1)),vox(3)*(1:size(skel,3)));

fig4 = figure(4);
col=[.7 .7 .8];
col2 = [1 0 0];
hiso = patch(isosurface(X,Y,Z,aortaMask - maskedVolume,0),'FaceColor',col,'EdgeColor','none');
%hiso = patch(isosurface(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isosurface(X,Y,Z,maskedVolume,0),'FaceColor',col2,'EdgeColor','none');
axis equal;axis off;
lighting gouraud;
isonormals(X,Y,Z,aortaMask,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
isonormals(X,Y,Z,maskedVolume,hiso2);
[x,y,z]=ind2sub([w,l,h],find(labeledskel>0));
plot3(vox(2)*y,vox(1)*x,vox(3)*z,'square','Markersize',4,'MarkerFaceColor','k','Color','k');
plot3(junctionCoords(2),junctionCoords(1),junctionCoords(3),'o','Markersize',8,'MarkerFaceColor','y','Color','k') 
set(gcf,'Color','white');
view(140,80)
title(sprintf('Ascending aorta volume %g mL',AAvolume/1000))
% Save the figure
saveas(fig4,[results_folder filesep() 'skel.png'])

% %% Plot the planes
% % Generate the patches for plotting the planes
% % Preallocate
% Xpatch = zeros(4,numberOfPlanes);
% Ypatch = zeros(4,numberOfPlanes);
% Zpatch = zeros(4,numberOfPlanes);
% for ii = 1:numberOfPlanes
%    Xpatch(1,ii) = planeStruct(ii).coords(1,1,1);
%    Xpatch(2,ii) = planeStruct(ii).coords(1,end,1);
%    Xpatch(3,ii) = planeStruct(ii).coords(end,end,1);
%    Xpatch(4,ii) = planeStruct(ii).coords(end,1,1);
%    
%    Ypatch(1,ii) = planeStruct(ii).coords(1,1,2);
%    Ypatch(2,ii) = planeStruct(ii).coords(1,end,2);
%    Ypatch(3,ii) = planeStruct(ii).coords(end,end,2);
%    Ypatch(4,ii) = planeStruct(ii).coords(end,1,2);
%    
%    Zpatch(1,ii) = planeStruct(ii).coords(1,1,3);
%    Zpatch(2,ii) = planeStruct(ii).coords(1,end,3);
%    Zpatch(3,ii) = planeStruct(ii).coords(end,end,3);
%    Zpatch(4,ii) = planeStruct(ii).coords(end,1,3);
% end
% figure();
% col=[.7 .7 .8];
% col2 = [1 0 0];
% hiso = patch(isosurface(X,Y,Z,aortaMask - maskedVolume,0),'FaceColor',col,'EdgeColor','none');
% hiso2 = patch(isosurface(X,Y,Z,maskedVolume,0),'FaceColor',col2,'EdgeColor','none');
% axis equal;axis off;
% lighting gouraud;
% isonormals(X,Y,Z,aortaMask,hiso);
% alpha(0.5);
% set(gca,'DataAspectRatio',[1 1 1])
% camlight;
% hold on;
% isonormals(X,Y,Z,maskedVolume,hiso2);
% patch(Ypatch,Xpatch,Zpatch,linspace(0,1,numberOfPlanes))
% set(gcf,'Color','white');
% view(140,80)

%% Plot the surfaces
for tt = tsystole
%for tt = 1:size(planeStruct(1).Vparallel,3)
    fig5 = figure(5);
    col=[.7 .7 .8];
    hiso = patch(isosurface(X,Y,Z,aortaMask,0),'FaceColor',col,'EdgeColor','none','FaceAlpha',0.2);
    axis equal;axis off;
    set(gca,'DataAspectRatio',[1 1 1])
    hold on;
    for ii = 1:size(planeStruct,2)
        if(size(planeStruct(ii).coords,1) > 1)
            surf(planeStruct(ii).coords(:,:,2),planeStruct(ii).coords(:,:,1),planeStruct(ii).coords(:,:,3),planeStruct(ii).Vparallel(:,:,tt),'EdgeColor','none');
        end
        %surf(planeStruct(ii).coords(:,:,2),planeStruct(ii).coords(:,:,1),planeStruct(ii).coords(:,:,3),planeStruct(ii).Vparallel(:,:,tt)./max(max(planeStruct(ii).Vparallel(:,:,tt))),'EdgeColor','none');
        %surf(planeStruct(ii).coords(:,:,2),planeStruct(ii).coords(:,:,1),planeStruct(ii).coords(:,:,3),planeStruct(ii).Vy(:,:,tt)./max(max(planeStruct(ii).Vy(:,:,tt))),'EdgeColor','none');
    end
    set(gcf,'Color','white');
    set(gca,'Ydir','reverse')
    colorbar
    %view(140,80)
    view(2)
    title('Planes with velocity at systole')
    % Save the figure
    saveas(fig5,[results_folder filesep() 'planes_t=' num2str(tt) '.png'])
end


%% Plot the major and minor axes from the planes
xvals = linspace(0,100,numberOfPlanes);
majory = 0.1*cat(1,planeStruct.MajorAxis); % 0.1 for mm --> cm
minory = 0.1*cat(1,planeStruct.MinorAxis); % 0.1 for mm --> cm
fig6 = figure(6);
plot(xvals,majory,xvals,minory)
legend('Major axis','Minor Axis')
xlabel('Percent distance along AA (%)')
ylabel('Length (cm)')
title('Major and Minor Axes On Each Placed Plane')
% Save the figure
saveas(fig6,[results_folder filesep() 'majorAndMinorAxis.png'])

fprintf('Code execution complete in %i seconds.\n\n',round(toc));
fprintf('================================================:)\n\n\n');

% % Plot the mask
% % Find the indicies of the 
% ind = find(maskedVolume == 1);
% [x,y,z] = ind2sub(size(labeledskel),ind);
% x = x.*vox(1);
% y = y.*vox(2);
% z = z.*vox(3);
% ind2 = find(aortaMask == 1);
% [x2,y2,z2] = ind2sub(size(labeledskel),ind2);
% x2 = x2.*vox(1);
% y2 = y2.*vox(2);
% z2 = z2.*vox(3);
% figure()
% scatter3(x,y,z)
% hold on
% scatter3(x2,y2,z2,'*r')
% axis('equal')

%     theta = linspace(0,2*pi,12).';
%     normalPlanes{ii,1} = bsxfun(@plus,R*(cos(theta)*T(1,:)+sin(theta)*T(2,:)),P0);
%     planePatchData{ii} = planePatchGenerator(interpolatedMidline,indicies(ii),interpolatedVelData,newGridVectors,timeStep);
%     disp([' - Plane ' num2str(ii) '/' num2str(length(indicies)) ' analyzed'])

% % Plot the normal vectors on the midline
% % Generate a meshgrid
% [X,Y,Z] = meshgrid(vox(2)*(1:size(skel,2)),vox(1)*(1:size(skel,1)),vox(3)*(1:size(skel,3)));
% 
% for ii = 1:size(planeStruct,2)
%     qq = planeStruct(ii).coords;
%     xq(ii) = qq(ceil(size(qq,1)/2),ceil(size(qq,2)/2),1);
%     yq(ii) = qq(ceil(size(qq,1)/2),ceil(size(qq,2)/2),2);
%     zq(ii) = qq(ceil(size(qq,1)/2),ceil(size(qq,2)/2),3);
%     uq(ii) = planeStruct(ii).Nhat(1);
%     vq(ii) = planeStruct(ii).Nhat(2);
%     wq(ii) = planeStruct(ii).Nhat(3);
% end
% 
% figure();
% col=[.7 .7 .8];
% col2 = [1 0 0];
% hiso = patch(isosurface(X,Y,Z,aortaMask - maskedVolume,0),'FaceColor',col,'EdgeColor','none');
% hiso2 = patch(isosurface(X,Y,Z,maskedVolume,0),'FaceColor',col2,'EdgeColor','none');
% axis equal;axis off;
% lighting gouraud;
% isonormals(X,Y,Z,aortaMask,hiso);
% alpha(0.5);
% set(gca,'DataAspectRatio',[1 1 1])
% camlight;
% hold on;
% isonormals(X,Y,Z,maskedVolume,hiso2);
% quiver3(yq,xq,zq,vq,uq,wq)
% set(gcf,'Color','white');
% view(140,80)

end
