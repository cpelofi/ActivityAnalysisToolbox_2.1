function [Chi,pVal,DAct,Bins,v1,v2] = relatedActivitiesTest(AllC1,AllC2,NBins)

% function [Chi,pVal,DAct,Bins,v1,v2] = relatedActivitiesTest(AllC1,AllC2,k)
% This function performes a goodness of fit test on the joint activity 
% distribution of two activities assessed on matching time frames 
% against the assumption that these are independent. 
%
% AllC* are matricies of binary point process, one column per time series
% marking which time frames contained events. These can be generated using
% activityCount.m: [~,AllC]=activityCount(Time,Series,FrameSize,HopSize,Thresh,option)
% For this test, each Series must be of the same length, aligned to the 
% same stimulus, and with FrameSize = HopSize.
% NBins is the prefered number of bins (per row and per columns) for the
% contigency table used for this test. NBins defaults to 3 if it isn't
% specified in the input arguments.
%
% The output Dact is a 3 dimensional array: 
% Dact(:,:,1) is the actual distribution for combinations of activity-levels
% Dact(:,:,2) is the estimated distribution if the series were independent.
% Bins is a similar kXkX2 array for the (actual, model) counts per bin used
% in the goodness of fit test. k is reduced to a min of 2 if any estimated
% bins contain less than 5 samples.
% Chi is the chi squared value from comparing the two sets of bins. 
% pVal is the likelyhood of the actual distribution being the result of two
% independent processes. 

% This function requires the standard MatLab Statistics Toolbox (chi2cdf.m)
% and the Activity toolbox function equiSplit.

% Finn Upham, April 6th, 2012
% Updated 2012/08/23
% Updated 2016/01/27, renamed, previously part of function jointChiSq

Chi = 0;
pVal = 0;
DAct = [];
Bins = [];
v1 = [];
v2 = [];

AC1 = sum(AllC1,2);
N1 = size(AllC1,2);

AC2 = sum(AllC2,2);
N2 = size(AllC2,2);

if nargin<3
    NBins = 3;
end

if size(AllC1,1) < 50
    NBins = 2;
end

L=size(AC1,1);
if size(AC2,1)~=L
    error('MatLab:jointChiSq',...
        'The activity time series are not of matching length')
end



% Calculate the joint distributions, actual and model.

% independent distribution, for independent activity
% test the two activities per time frame against
% a null hypothesis of independence. These activities must not be
% exclusive.

% joint distribution of activities on time frames for all activity-
% levels
A=zeros(N1+1,N2+1);
for i=0:N1
    A(i+1,:)=hist(AC2(AC1==i),0:N2);
end

% hypothetical independent distribution of activities for k ranges 
% of activity levels.
% Activity 1
[NAct1]=hist(AC1,0:N1);
[v1,bins1] = equiSplit(NAct1,NBins,NBins*6);
% Activity 2
[NAct2]=hist(AC2,0:N2);
[v2,bins2] = equiSplit(NAct2,NBins,NBins*6);

% expected distribution over all joint activity-levels
B = NAct1'*NAct2/L;

% dof = 0; % degrees of freedom adjustment

% actual distribution
A=zeros(N1+1,N2+1);
for j=0:N1
    A(j+1,:)=hist(AC2(AC1==j),0:N2); %make sure this isn't backwards
end


% In case the binnings of each separate dimension did not yeild sufficient
% expected samples in each table entry, we reduce the number of rows or
% columns accordingly.

%check out the models table to see if the binnings in both dimensions
%resulted in appropriate distribution of samples in all entries of the
%table
tableB = zeros(length(bins1),length(bins2));

if isempty(tableB)
    Chi = 0;
    pVal = 1;
    DAct = [];
    Bins = cell(0);
    v1 = [];
    v2 = [];
    return;
end
    
for i = 1:length(bins1)
    for j = 1:length(bins2)
        tableB(i,j) = sum(sum(B(v1{i},v2{j})));
    end
end

%check the quality of the model distribution, shrink the table if bad.
while min(min(tableB))<5
    if max(size(tableB))>2
        if length(bins2)<length(bins1)
            [v1,bins1] = equiSplit(NAct1,NBins-1,16);
        elseif length(bins2)>length(bins1)
            [v2,bins2] = equiSplit(NAct2,NBins-1,16);
        elseif  min(bins1)<min(bins2)
            [v1,bins1] = equiSplit(NAct1,NBins-1,16);
            NBins = NBins-1;
        else
            [v2,bins2] = equiSplit(NAct2,NBins-1,16);
            NBins = NBins-1;
        end

        tableB = zeros(length(bins1),length(bins2));
        for i = 1:length(bins1)
            for j = 1:length(bins2)
                tableB(i,j) = sum(sum(B(v1{i},v2{j})));
            end
        end
    else
        %fprintf('Bad distributions, can not test these activity levels\n')
        pVal = 1; %Nan
        Chi = 0; %Nan
        return;
    end
end


% table for actual
tableA = zeros(size(tableB));
for i = 1:length(bins1)
    for j = 1:length(bins2)
        tableA(i,j) = sum(sum(A(v1{i},v2{j})));
    end
end


DAct(:,:,1)=B;
DAct(:,:,2)=A;

% Finally calculate the Goodness of fit test

Bins=zeros(length(bins1),length(bins2),2);
Bins(:,:,1)=tableA;
Bins(:,:,2)=tableB;

Chi=sum(sum(((Bins(:,:,2)-Bins(:,:,1)).^2)./Bins(:,:,2)));

pVal=1-chi2cdf(Chi,length(bins1)+length(bins2)-2);



 