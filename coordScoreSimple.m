function C = coordScoreSimple(Time,Series,winSize,Thresh,option,Nbins)

%  C = coordScoreSimple(Time,Series,winSize,Thresh,option,Nbins)
%
% This function calculates the composite coordination score for a Series 
% for a specified activity event type over non-overlapping time frames. 
% The type of activity is defined by option and Thresh as inputs to
% activityCount.m.
%
% This calculates the goodness of fit p value for winSize segmentations of
% the time series into winsize timeframes by shifting the starting point by 
% a sample per iteration. This shifting captures concentrations of activity
% levels which might otherwise be split by some segmentations.
% these pvalues are aggregated into a score: the mean -log_10 of these values.

% Nbins is the optional input defining the maximum number of divisions of the
% activity-level distributions used for the goodness of fit test. If Nbins
% isn't specified, the default is 4.

% C is the coordination score calculated on the parametric test pvalues

% This function uses activityCount and simpleActivityTest from the Activity
% Analysis Toolbox

% Finn Upham
% reviewed 2014/07/06
% reviewed and renamed 2016 01 27, previously monoCoord

if nargin < 6
   Nbins = 4;
end

CS = zeros(winSize,1);

for i=1:winSize
    AC=activityCount(Time(i:end),Series(i:end,:),winSize,winSize,Thresh,option);
    [~,CS(i)]=simpleActivityTest(AC,size(Series,2),Nbins);
end
 
C=mean(-log10(CS+10^(-16)));