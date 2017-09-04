function [veloMIP,vmaxcoords,tsystole] = velocityMIP(vel,mask)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% % Erode segmentation to offset segmentation errors
se = strel(ones(3,3,3));
erodedMask = imerode(mask,se);

%erodedMask = mask;

% Calculate the velocity magnitude at each voxel
velocityMagnitude = zeros(size(vel,1),size(vel,2),size(vel,3),size(vel,5));
for tt = 1:size(vel,5)
    % Apply the mask
    vel(:,:,:,1,tt) = vel(:,:,:,1,tt) .* erodedMask;
    vel(:,:,:,2,tt) = vel(:,:,:,2,tt) .* erodedMask;
    vel(:,:,:,3,tt) = vel(:,:,:,3,tt) .* erodedMask;
    
    velocityMagnitude(:,:,:,tt) = (vel(:,:,:,1,tt).^2 + vel(:,:,:,2,tt).^2 + vel(:,:,:,3,tt).^2).^0.5;
end

% Find systole
averageVelocity = zeros(1,size(vel,5));
for tt = 1:size(vel,5)
    averageVelocity(tt) = sum(sum(sum(velocityMagnitude(:,:,:,tt)))) / sum(sum(sum(velocityMagnitude(:,:,:,tt) ~= 0)));
end

[~,idx] = max(averageVelocity);
if idx == 1
    systoleIdx = [1 2 3];
elseif idx == length(averageVelocity)
    systoleIdx = [idx-2 idx-1 idx];
else 
    systoleIdx = [idx-1 idx idx+1];
end
tsystole = idx;
clear idx;

    
% Grab velocity values from the eroded segmentation at systole +/- 1 time
% point. Get only the maximum over the time point
systoleVelocities = max(velocityMagnitude(:,:,:,systoleIdx),[],4);
orderedVelocities = nonzeros(systoleVelocities);
% Order the values from lowest to highest
orderedVelocities = sort(orderedVelocities);

% Take the diff of the velocity values
D = diff(orderedVelocities);
% T = C * mean (D), C = 10; from Rose, et al.
C = 15;
T = C * mean(D);

% Find the index of the first value that is over T
idx = find(D > T);
% A few noisy low values cause problems, look only in the second half of D
idx = min(idx(idx > length(D) / 2));
% Eliminate all velocities over the value of T
filteredVelocities = orderedVelocities(1:idx-1); 
% Find the maximum velocity value idx
vmaxidx = find(velocityMagnitude == filteredVelocities(end));
% Convert to coordinates
[vmaxx,vmaxy,vmaxz,vmaxt] = ind2sub(size(velocityMagnitude),vmaxidx);
vmaxcoords = [vmaxx vmaxy vmaxz];
% What to do if there are two maxima?

systoleVelocities(systoleVelocities == 0) = NaN;

% Generate the MIPs
MIP1 = squeeze(max(systoleVelocities,[],1));
MIP2 = squeeze(max(systoleVelocities,[],2));
MIP3 = max(systoleVelocities,[],3);

veloMIP.MIPz = MIP3;
veloMIP.MIPzcoords = [vmaxx vmaxy];
veloMIP.MIPy = MIP2;
veloMIP.MIPycoords = [vmaxx vmaxz];
veloMIP.MIPx = MIP1;
veloMIP.MIPxcoords = [vmaxy vmaxz];
end

