function [searchData] = calculatePaths(node, searchData, currentPath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Check if the next node has branches
currentNode = currentPath(end);
if(node(currentNode).ep == 0 || length(currentPath) == 1)
    % Not an endpoint
    for jj = 1:length(node(currentNode).conn)
        % Check if the next node has already been traversed in this path
        if(~any(currentPath == node(currentNode).conn(jj)))
            newCurrentPath = [currentPath node(currentNode).conn(jj)];
            searchData = calculatePaths(node, searchData, newCurrentPath);
        end
    end
elseif(node(currentNode).ep == 1)
    % Is an endpoint
    %completePath = [currentPath node(currentNode).conn(1)];
    % Write to the structure
    searchData(end+1).path = currentPath;
else
    error('Unexpected value: node.ep should be 0 if not an endpoint or 1 if an endpoint');
end
end


