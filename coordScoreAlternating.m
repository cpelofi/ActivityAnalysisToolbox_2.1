function CS = coordScoreAlternating(Time,Series,winSize,ThreshA,optionA,ThreshB,optionB,k)

% function CS = coordScoreAlternating(Time,Series,winSize,ThreshA,optionA,ThreshB,optionB,k)
% this function calculates the parametric coordination score for two types 
% of activity in a single collection of time series Series, assessed on 
% single activities over non-overlapping time frames.

% It calculates the p value for the contingency tables (independence) by shifting
% the cutting of time frames by i=1:winsize samples before assessing the activity
% and reports the mean -log_10 of these pvalues.
% CS is the coordination score calculated on the parametric test pvalues

% this function uses others in the Activity Analysis Toolbox:
% actionCount and alternatingActivityTest

% Finn Upham
% reviewed 2014/07/06
% reviewed 2016 - 01 - 27, renamed from altCoord


if nargin < 9
    k = 3;
end

CS = zeros(winSize,1);

for i=1:winSize

    [~,AllCa]=activityCount(Time(i:end),Series(i:end,:),winSize,...
        winSize,ThreshA,optionA);
    [~,AllCb]=activityCount(Time(i:end),Series(i:end,:),winSize,...
        winSize,ThreshB,optionB);
    
    [~,CS(i)]=alternatingActivitiesTest(AllCa,AllCb,k);
end

CS=mean(-log10(CS+10^(-16)));