function [Vinterp] = interpolateVelocity(X,Y,Z,V,Q,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% To do:
% - Implement a check to see if the interpolation is out of bounds. What to
% do in this case? NaN? Would need to implement error handling in the
% subsequent code

%% Check inputs
if nargin < 6
    method = 'linear';
end

% If the method is not one of the accepted methods, default to 'linear'
if nargin == 6
    method = lower(method);
    if sum(strcmp(method,{'nearest', 'cubic', 'spline'})) == 0
        method = 'linear';
    end
end

%% Convert the query points into x,y,z
xvec = zeros(size(Q,1)*size(Q,2),1);
yvec = xvec;
zvec = xvec;

index = 1;
for jj = 1:size(Q,2)
    for ii = 1:size(Q,1)
        xvec(index) = Q(ii,jj,1);
        yvec(index) = Q(ii,jj,2);
        zvec(index) = Q(ii,jj,3);
        index = index + 1;
    end
end

%% Make the interpolant
P = [2 1 3];
X = permute(X, P);
Y = permute(Y, P);
Z = permute(Z, P);
V = permute(V, P);
F = griddedInterpolant(X,Y,Z,V,method);

%% Interpolate
Vinterp = zeros(size(Q,1),size(Q,2));
for ii = 1:length(xvec)
    % Wasn't working for F(xvec(ii),xvec(ii),zvec(ii))
    Vinterp(ii) = F(yvec(ii),xvec(ii),zvec(ii));
end
end

