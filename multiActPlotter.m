function [acts,dT]=multiActPlotter(Time,Data,winSize,hopSize,option,thresholds)

%function [acts,dT,h]=multiActPlotter(Time,Data,winSize,hopSize,option,intervals)
% generates activity levels for multiple thresholds for events in Data
% specified by option. If no outputs are specified, the function also 
% plots these in an area styles figure.

% Input variables correspond to those of activityCount.m

% Finn Upham 2016 01 23

thresholds = sort(thresholds,'descend');
ACs = [];

% actionCount can't handle nan's, so out they go.
ind = sum(Data,2);
Data = Data(~isnan(ind),:);
Time = Time(~isnan(ind));

for i = 1:length(thresholds)
    [AC,~,dT]=activityCount(Time,Data,winSize,hopSize,thresholds(i),option);
    ACs = [ACs AC];
end

acts = [ACs(:,1) diff(ACs')'];

% the area plot 
if nargout < 2
    figure
    acts = area(dT,acts);
    axis([Time(1) Time(end) 0 1])
end

