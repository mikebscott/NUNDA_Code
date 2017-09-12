function [max_vel, xcoord, ycoord, zcoord] = velMIPfilter(velmag_all)

% Find the size of the velocity data
velSize = size(velmag_all);
% if length(velSize) == 3
%     % Do nothing
% elseif length(velSize) == 4
%     velSize = velSize(1:3);
% else
%     error('velmag_all is the wrong size')
% end

velmag_all(velmag_all == 0) = NaN;

%Squeeze velmag into vector-> Sort Vector(S_vector)--> Remove NAN values.
velVector = velmag_all(:);

% Find the maximum value of the velocity
max_val = max(velVector);
% Find the idx of the vmax
vmaxidx = find(velmag_all == max_val);
% Convert to subscripts
[xcoord, ycoord, zcoord] = ind2sub(velSize,vmaxidx);

%Filter Threshold Scale 10 for velcotiy 100 for KE
mlt = 10;

filter10_FLAG = [];

% Old for loop
[S_vector, IND] = sort(velVector);
S_vector(isnan(S_vector)) = [];
S_vectorL = length(S_vector);

%If S_vector Length is greater than 3000 take top 0.5% values. If not take top 15.
if S_vectorL < 3000
    toptierL = 15;
else
    toptierL = round(0.005*S_vectorL);
end

%Get median of top values for velocity.
toptier = S_vector(S_vectorL-toptierL+1: S_vectorL);
if mod(toptierL,2) == 0
    median_toptier = toptier((toptierL/2+1));
else
    median_toptier = median(toptier);
end

%% MeanFilter: Takes derivative (d1) of top 1% values. Starting at the bottom of d1, if any value exceeds 10*mean of d1 that value is taken as the max velocity.
toptierL_2 = round(0.01*S_vectorL);
toptier_2 = S_vector(S_vectorL-toptierL_2+1: S_vectorL);
d1 = diff(toptier_2);

for ii = 1:length(d1)
    if d1(ii) >= mlt*mean(d1)
        max_vel=toptier_2(ii);
        filterIND=(S_vectorL-toptierL_2+ii);
        filter10_FLAG=1;    %Flags that filter has excluded voxel(s)
        break; % Break if the derivative value exceeds the threshold!
    end
end
if exist('max_vel', 'var')==0
    max_vel = max_val;
end


%If the max velocity given from MeanFilter is greater than 1.2*median
%median is used as the new max velocity.

if max_vel > 1.2*median_toptier
    max_vel = median_toptier;
    filterIND=find(S_vector == median_toptier, 1);
    [xcoord, ycoord, zcoord] = ind2sub([velSize(1) velSize(2) velSize(3)], IND(filterIND));
elseif isempty(filter10_FLAG)==0
    [xcoord, ycoord, zcoord] = ind2sub([velSize(1) velSize(2) velSize(3)], IND(filterIND));
end

if isempty(xcoord) == 1
    % Use the maximum velocity
    idx = find(velmag_all == max_vel);
    [xcoord, ycoord, zcoord] = ind2sub(velSize,idx);
end

end