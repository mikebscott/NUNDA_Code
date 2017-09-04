function searchData = calculatePathlengths(searchData,A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% For each of the different paths found
for ii = 1:size(searchData,2)
    % Initialize the distance
    distance = 0;
    % Calculate the distance along the path
    for jj = 2:length(searchData(ii).path)
       distance = distance + A(searchData(ii).path(jj-1),searchData(ii).path(jj)); 
    end
    % Save the distance in the searchData struct
    searchData(ii).distance = distance;
end

end