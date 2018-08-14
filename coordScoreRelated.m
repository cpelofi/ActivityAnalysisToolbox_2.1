function CS = coordScoreRelated(Time,SeriesA,SeriesB,winSize,Thresh1,option1,Thresh2,option2,Nbins)

% function CS = coordScoreRelated(Time,SeriesA,SeriesB,winSize,Thresh1,option1,Thresh2,option2,Nbins)
% this function calculates the parametric coordination score for two types 
% of activity in two synchronously sampled collections of time series 
% SeriesA and SeriesB, assessed on single activities (Thresh1, option1, 
% Thresh2, option2) over non-overlapping time frames.
% It calculates the p value for the contingency tables (independence) by shifting
% the cutting of time frames by i=1:winsize samples before assessing the activity
% and reports the mean -log_10 of these pvalues.
% cSa is the coordination score calculated on the parametric test pvalues

% this function uses actionCount and jointChiSq

% Finn Upham
% reviewed 2014/07/06
% reviewed 2016/03/22

CS = zeros(winSize,1);

if nargin < 7
    Nbins = 3;
    option2 = option1;
    Thresh2 = Thresh1;
elseif nargin < 8
    Nbins = Thresh2;
    option2 = option1;
    Thresh2 = Thresh1;
elseif nargin < 9
    Nbins = 3;
end

% per slicing, calculate the activity
for i=1:winSize
    [~,AllCa]=activityCount(Time(i:end),SeriesA(i:end,:),winSize,...
        winSize,Thresh1,option1);
    [~,AllCb]=activityCount(Time(i:end),SeriesB(i:end,:),winSize,...
        winSize,Thresh2,option2);
    
    [~,CS(i)]=relatedActivitiesTest(AllCa,AllCb,Nbins);
end

CS=mean(-log10(CS+10^(-16)));