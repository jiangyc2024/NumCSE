
% Problem 5g)
% Load the data
load fotodaten.dat;

% Set the initial values of the counters to 0
StarCnt=0;
XmasCnt=0;

% Compute the original points of the star and of the xmas tree
StarPts = star;
XmasPts = Xmastree;
close all;

for i = 1 : 20
    % Load the i-th collection of 15 points denoting a figure
    PtsFromData = fotodaten(2*i-1:2*i, :);
    XmasCn
    % Compute both errors
    [~, err1] = complinmap(StarPts, PtsFromData);
    [~, err2] = complinmap(XmasPts, PtsFromData);
    
    % Error comparison yields the results
    if err1 > err2 
        XmasCnt = XmasCnt+1;
    else
        StarCnt = StarCnt+1;
    end
end