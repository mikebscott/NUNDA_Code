function [planeStruct] = planeAnalysis(planeStruct,tr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nPlanes = size(planeStruct,2);
tsteps = size(planeStruct(1).Vx,3);
if(tr > 1 && tr < 200)
    tr = tr / 1000; % s from ms
elseif(tr >= 200)
    % Assume tr is the nominal interval
    nominalInterval = tr/1000;
end

% Calculate the resolution of the plane (in mm)
for ii = 1:size(planeStruct,2)
    if(size(planeStruct(ii).coords,1) > 1)
        pix = sqrt(sum((planeStruct(ii).coords(1,1,:)-planeStruct(ii).coords(1,2,:)).^2));
        A = pix^2; % in mm^2
        break;
    end
end


% Calculate peak flow
for ii = 1:nPlanes
    % Preallocate storage variables
    netflow = zeros(1,tsteps);
    retrogradeflow = zeros(1,tsteps);
    forwardflow = zeros(1,tsteps);
    for tt = 1:tsteps
        pixelflow = planeStruct(ii).Vparallel(:,:,tt).*A;
        pixelflow(isnan(pixelflow)) = 0;
        netflow(tt) = sum(sum(pixelflow));
        retrogradeflow(tt) = sum(sum(pixelflow(pixelflow < 0)));
        forwardflow(tt) = sum(sum(pixelflow(pixelflow > 0)));
        %retrogradefraction = abs(retrogradeflow) / forwardflow;
    end
    
    % Store the outputs
    planeStruct(ii).netflow = netflow;
    planeStruct(ii).forwardflow = forwardflow;
    planeStruct(ii).retrogradeflow = retrogradeflow;
    % Cardiac output should be forward flow?
    if(exist('nominalInterval','var'))
        planeStruct(ii).cardiacoutput = sum(planeStruct(ii).forwardflow)*60/1000/nominalInterval/100; % L/min       
    else
        planeStruct(ii).cardiacoutput = sum(planeStruct(ii).forwardflow)*60/1000/(tsteps*tr)/100; % L/min
    end
    % Stroke volume, define as volume ejected during systole
    % Find the idx of the first negative net flow (assume start diastole)
    idx = find(planeStruct(ii).netflow < 0,1);
    if(idx == 1)
        temp = find(planeStruct(ii).netflow < 0,2);
        if(length(temp) == 2)
            idx = temp(2);
        else
            idx = 1;
        end
        clear temp;
    end
    planeStruct(ii).strokevolume = sum(netflow(1:idx-1))*tr;
    
    planeStruct(ii).peakvelocity = max(max(max(planeStruct(ii).Vparallel)));
    % Find the index corresponding to the peak velocity
    idx = find(planeStruct(ii).Vparallel == planeStruct(ii).peakvelocity);
    [x, y, z] = ind2sub(size(planeStruct(ii).Vparallel),idx);
    planeStruct(ii).peakvcoords = [x y z];
    % Define retrograde fraction as retrograde flow / net flow
    planeStruct(ii).retrogradefraction = abs(sum(planeStruct(ii).retrogradeflow)) / sum(planeStruct(ii).netflow);
end





% % plot
% figure()
% hold on
% for ii = 1:nPlanes
%     plot(1:tsteps,planeStruct(ii).netflow,'-k')
%     plot(1:tsteps,planeStruct(ii).forwardflow,'-b')
%     plot(1:tsteps,planeStruct(ii).retrogradeflow,'-r')
% end
end