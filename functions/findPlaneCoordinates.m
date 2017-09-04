function [planeCoordinates] = findPlaneCoordinates(P0, Nhat, maxWidth, spacing)
% Function

if(nargin < 4)
   spacing = 1; % Default 1 mm 
end

% Find the number of elements needed to cover the space
numberOfVoxels = ceil(maxWidth / spacing);
if mod(numberOfVoxels,2) == 0
    % If even, make it odd
    numberOfVoxels = numberOfVoxels + 1;
end

% Generate a spacing vector
spacingVector = (-1*floor(numberOfVoxels/2)*spacing:spacing:floor(numberOfVoxels/2)*spacing);

% Allocate the matrix to hold the coordinates
%   x*y*3 matrix
%   dim 1: x coordinate
%   dim 2: y coordinate
%   dim 3: z coordinate
planeCoordinates = zeros(numberOfVoxels,numberOfVoxels,3);
planeCoordinates(:,:,1) = repmat(spacingVector,numberOfVoxels,1);
planeCoordinates(:,:,2) = repmat(spacingVector,numberOfVoxels,1)';

% Find the axis of rotation
% Note that the plane being rotates is perpendicular to the z axis, so
% we're looking for the rotational axis between the normal vector input and
% the z axis.
Nz = [0 0 1];
R = RotationFromTwoVectors(Nz,Nhat);

% Check if the vectors are already aligned

for ii = 1:numberOfVoxels
    for jj = 1:numberOfVoxels
        % Apply the rotation and translate by P0
        coords = [planeCoordinates(ii,jj,1) planeCoordinates(ii,jj,2) planeCoordinates(ii,jj,3)]; 
        newCoords = R*coords';
        planeCoordinates(ii,jj,1) = newCoords(1) + P0(1);
        planeCoordinates(ii,jj,2) = newCoords(2) + P0(2);
        planeCoordinates(ii,jj,3) = newCoords(3) + P0(3);
    end
end

% Check for any 0 or negative coordinates, replace with NaN!
end

function R=RotationFromTwoVectors(A, B)
% Aligns vector A with vector B
v = cross(A,B);
ssc = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + ssc + ssc^2*(1-dot(A,B))/(norm(v))^2;
end