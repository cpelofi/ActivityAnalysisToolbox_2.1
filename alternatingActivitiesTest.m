function [Chi,p,DAct,Bins,v1,v2] = alternatingActivitiesTest(AllC1,AllC2,k)

% function [Chi,p,DAct,Bins,v1,v2] = alternatingActivitiesTest(AllC1,AllC2,k)
% This function performes a goodness of fit test on the activity 
% distribution of two activities evaluated on the same collection of time 
% series, such that these two activities can not both occure in the same
% time frame.
% They are tested against a null hypothesis of independence between their
% combined activity levels and the ratio of between these two activities
% per frame.

% For example, increases and decreases in continuous ratings, as defined 
% in actionCount, are mutually exclusive within a time frame. If the two
% forms of activity are alternating, moments with high activity-levels will
% be predominently of one type. If, however, they are independent of each
% other, moments of high activity will tend to be evenly split. The null
% hypothesis is derived from the actual activity-levels distributions 
% for both types of activity.

% The output Dact is a 3D array with the first matrix being the actual
% joint distribution for all numbers of participants, and the second
% being the estimated distribution if the series were independent.
% Bins is 3X3X2 array as well (actual, model) of the bins selected for the
% goodness of fit test. If the bins hold too few samples (< 5), the bins 
% are recalculated to generate 6 bins rather than 9 to describe the
% distribution.

% Chi is the chi squared value from comparing the two sets of bins. 

% p is the likelyhood of the actual distribution being the result of
% independent processes. This requires the Statistics Toolbox. Also depends
% on function equiSplit.

% Finn Upham, April 6th, 2012
% Updated 2012/08/23 Alternating solved!
% Updated 2016/01/27, renamed, previously part of function jointChiSq

% error checking
% matching size
L=size(AllC1);
if size(AllC2,1)~=L(1)
    error('MatLab:alternatingActivityTest',...
        'The activity event series do not report the same number of time frames. The two collections of activity point processes must report events in matching time frames.')
end

if size(AllC2,2)~=L(2)
    error('MatLab:alternatingActivityTest',...
        'The activity event series do not report the same number of time series. The two collections of activity point processes must report events occuring in the same time series.')
end

% insufficient number of frames for a 3 by 3 table of activity levels
if L(1) < 50
   error('MatLab:alternatingActivityTest',...
        'The activity event series do not report a sufficent number of time frames to apply this test.')
end

% check if activity events are exclusive
if (sum(AllC1(logical(AllC2)))+sum(AllC2(logical(AllC1))))>0
    error('MatLab:alternatingActivityTest',...
        'These two Activities are not exclusive, thus this test can not be applied.')
end


% preping variables
Chi = 0;
p = 0;
DAct = [];
Bins = [];
v1 = [];
v2 = [];

if nargin<3
    k = 3;
end

% Calculate the joint distributions, actual and model.
% this is why I hate combinatorics
% new, shinier and hopefully more reasonable solution 
% (August 22, 2012)

AllC = sign(AllC1+AllC2);
N = size(AllC,2); %number of series in collection
L = size(AllC,1); %number of samples per series

AC = sum(AllC,2); %total activity count time series
p = sum(AC)/(N*L); %rate of any activity
NAct = hist(AC,0:N); %Actual distribution of activity levels (all activity)
ACa = sum(AllC1,2); % activity a count time series
pa = sum(ACa)/(N*L); % rate of activity a (event)
ACb = sum(AllC2,2); % activity b count time series
pb = sum(ACb)/(N*L); % (activity rate) estimate of likelyhood of event b
pab = pb/p; % if there is an event, the likelyhood of it being b

Model = zeros(N+1);
altModel = zeros(N+1);

for s = 0:N %activity sum
  set = zeros(s+1,1);
  if s < 50 % because cald
    for r = 0:s
      set(r+1) = nchoosek(s,r)*((1-pab)^(s-r))*(pab^r);
    end
  else
    for r=0:s
      set(r+1)=exp(-s*pab)*((s*pab)^r)/factorial(r);
    end 
  end

  set = NAct(s+1)*set/sum(set);
  altModel(1:s+1,s+1) = set;
  for r = 0:s
    Model(r+1,s+1-r) = set(r+1);
  end
end

B = Model;
A = zeros(N+1);
for j=0:N
    A(j+1,:)=hist(ACb(ACa==j),0:N); %make sure this isn't backwards
end
 
% mapping from s+1 values evenly to N+1 spots
mapper=cell(N+1,1);
mapper{1} = round(0.5*(N))+1;
for s = 1:N
    mapper{s+1} = round(N*(0:s)/s)+1;
end
 
% switch A and B to this form
Aalt = zeros(N+1);
Balt = Aalt;

for s = 0:N
  v = hist(ACb(AC==s),0:s);
  if ~isempty(v)
      Aalt(mapper{s+1},s+1)= v;
  end
  Balt(mapper{s+1},s+1)=altModel(1:s+1,s+1);
end

A = Aalt;
B = Balt;
 
% Now calculate the binning of the distributions

NAct1 = sum(B,2);
NAct2 = sum(B,1);

[v1,bins1] = equiSplit(NAct1,k,k*6);
[v2,bins2] = equiSplit(NAct2,k,k*6);
% dof = -1; % one (or two) estimated variables in alt distribution


% In case the binnings of each separate dimension did not yeild sufficient
% expected samples in each table entry, we reduce the number of rows or
% columns accordingly.

%check out the models table to see if the binnings in both dimensions
%resulted in appropriate distribution of samples in all entries of the
%table
tableB = zeros(length(bins1),length(bins2));

if isempty(tableB)
    C = 0;
    p = 1;
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
    if max(size(tableB))>3 % alternating needs a middle separate from the edges.
        if length(bins2)<length(bins1)
            [v1,bins1] = equiSplit(NAct1,k-1,16);
        elseif length(bins2)>length(bins1)
            [v2,bins2] = equiSplit(NAct2,k-1,16);
        elseif  min(bins1)<min(bins2)
            [v1,bins1] = equiSplit(NAct1,k-1,16);
            k = k-1;
        else
            [v2,bins2] = equiSplit(NAct2,k-1,16);
            k = k-1;
        end

        tableB = zeros(length(bins1),length(bins2));
        for i = 1:length(bins1)
            for j = 1:length(bins2)
                tableB(i,j) = sum(sum(B(v1{i},v2{j})));
            end
        end
    else
        %fprintf('Bad distributions, can not test these activity levels\n')
        p = 1; %Nan
        C = 0; %Nan
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

p=1-chi2cdf(Chi,length(bins1)+length(bins2)-2);



 