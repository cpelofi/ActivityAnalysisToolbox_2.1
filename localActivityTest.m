function [pVal,Coinc,CoincRank,CoincSurprise,AlternativeCoincs,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(AllC,FrameSize,ShuffleRange,Iter,option)

% function [pVal,Coinc,CoincRank,CoincSurprise,AlternativeCoincs,AltP,altpVal,NPCscore,altNPCscore]
% = localActivityTest(AllC,FrameSize,ShuffleRange,Iter,option)
%
% Version 3.0
%
% localActivity assess the likelyhood of the coincidence of events in a
% collection of time synchronised point processes in windows of size FrameSize,
% by non-parametrically estimating the distribution of coincidences by
% repeatedly randomly shifting the alignment of the series uniformly over
% time span tr.
%
% Inputs
% AllC: a collection of point process time series {0,1}, one series per column,
%       sampled at an constant sample rate and mutually aligned
%       Output of actionCount
% FrameSize: the window of coincidence, in samples of AllC
% ShuffleRange: the window of random shifting, in samples of AllC
% N: the interations of random shifting to inform the expecte distribution
% option: either 'Loop' or 'Reflect' this specifies the method of
%   approximating the expecte distribution for the ends of the series. Loop,
%   the default, pulls AllC from the opposite end of the series rather than
%   zero pad. Reflect fills AllC with the time reversed point process.
%
% Outputs
% pVal: the nonparametric p value of the coordination for the experimental 
%       data (aka rank amidst the Iter number of alternatives according to
%       the Eudlicidian distance metric from the average cumulative
%       distribution of activity-levels for the alternative alignments.)
% Coinc: a time series reporting the number of time series active in each
%       time frame, essentially the sum of AllC across columns.
%       The timestamps of the frames of Coinc align with those of AllC, as 
%       reported by the function activityCount.m
% CoincRank: the p value or percentile of the actual Coinc count per time
%       frame in the estimated distribution.
% CoincSurprise: the p-values translated into surprise values (-log10((1-p)/p))
% AlternativeCoincs: the activity-levels counted for each of the Iter time 
%       shuffling at each timeframe.
% AltP: Reported ranks of activity level distribution to the alternative 
%       distributions.
% altpVal: pValue of experimental data calculated from ranks distributions,
%       rather than activity level distributions of alternatives.
% NPCscore: The rank-based C Score for this nonparametric test of the
%       distribution of activity-levels. (C of pVal)
% altNPCscore: The rank-based C Score for this nonparametric test of the
%       distribution of local rank of activity-levels (a.k.a., surprise
%       distribution.) (C of altpVal)
% All of these outputs have the same number of rows as AllC.
% this method of estimating local surprise is derived from:
% Pipa, G., Wheeler, D., Singer, W., and Nikoli ?c, D. (2008). Neuroxidence:
%   reliable and efficient analysis of an excess or deficiency of joint-
%   spike events. Journal of computational neuroscience, 25(1):64?88.

% Finn Upham 2016 01 20
% Finn Upham 2017 05 30
% Finn Upham 2018 06 06 addition of Extend and Reflect options
% Finn Upham 2018 07 10 addition of altNPCscore
 
if nargin==4
    option='Loop';
end

L = size(AllC);
if round(ShuffleRange/2)~= ShuffleRange/2
    ShuffleRange = ShuffleRange+1;
end

% building alt AllC
AlternativeCoincs = zeros(L(1),Iter);

for j = 1:Iter
    shifts = round((rand(1,L(2)))*ShuffleRange);
    % extend the number of rows to include the tr shift range
    SAllC = zeros(L(1)+ShuffleRange,L(2));
    switch option
        case 'Ratio'  % calculate the incidence ratio on the smaller number of responses
            for i = 1:L(2)
                SAllC(shifts(i)+(1:L(1)),i) = AllC(:,i);
            end
            SAllC = SAllC((ShuffleRange/2)+1:end-(ShuffleRange/2),:);
            s = hist(shifts,1:ShuffleRange);
            s = cumsum(s);
            Z = ones(length(SAllC),1)*L(2);
            Z(1:ShuffleRange/2) = s((1+(ShuffleRange/2)):end);
            Z((1+end-ShuffleRange/2):end) = s(1:ShuffleRange/2);
            
            AC = sum(sign(conv2(ones(FrameSize,1),1,SAllC,'same')),2);
            AlternativeCoincs(:,j) = round(L(2)*AC./Z);
        case 'Extend' % extend from distribution at point of all coincidence
            for i = 1:L(2)
                SAllC(shifts(i)+(1:L(1)),i) = AllC(:,i);
            end
            SAllC = SAllC(ShuffleRange+1:end-(ShuffleRange),:);
            SC = sum(sign(conv2(ones(FrameSize,1),1,SAllC,'same')),2);
            SC = [repmat(SC(1,:),ShuffleRange/2,1);
                 SC;
                 repmat(SC(end,:),ShuffleRange/2,1)];
            AlternativeCoincs(:,j) = SC;
        case 'Reflect' % extend from distribution at point of all coincidence
            for i = 1:L(2)
                SAllC(shifts(i)+(1:L(1)),i) = AllC(:,i);
                SAllC(1:shifts(i),i) = AllC(shifts(i):-1:1,i);
                SAllC(end-(ShuffleRange-shifts(i))+1:end,i) = AllC(end:-1:end-(ShuffleRange-shifts(i))+1,i);
            end
            SAllC = SAllC((ShuffleRange/2)+1:end-(ShuffleRange/2),:);
            % count coincs over all ;
            AlternativeCoincs(:,j) = sum(sign(conv2(ones(FrameSize,1),1,SAllC,'same')),2);
        otherwise
            % case 'Loop'
            % loop uses the opposite end of the each series instead of zero padding
            % to fill the gaps left by offsetting the time series in each
            % iteration.
            
            for i = 1:L(2)
                d = AllC(:,i);
                D = SAllC(:,i);
                %[shifts(i) size(d) size(D)]
                % pad the begining and end of the shifted AllC with the
                % opposite end of the series
                % [i shifts(i) ShuffleRange size(d) size(D)]
                D(shifts(i)+(1:L(1))) = d;
                if shifts(i) < L(1)
                    D(1:shifts(i)) = d(end-shifts(i)+1:end);
                end
                if shifts(i) > 0
                    D(shifts(i)+L(1)+1:end) = d(1:ShuffleRange-shifts(i));
                end
                SAllC(:,i) = D;
                % [i shifts(i) tr length(d) length(d)-((tr/2)+shifts(i)-1) ((tr/2)-shifts(i))]
                %SAllC(shifts(i)+(1:L(1)),i) = [d(end-((ShuffleRange/2)+shifts(i)-1):end,1); AllC(:,i);d(1:((ShuffleRange/2)-shifts(i)),1)];
            end
            SAllC = SAllC((ShuffleRange/2)+1:end-(ShuffleRange/2),:);
            % count coincs over all ;
            AlternativeCoincs(:,j) = sum(sign(conv2(ones(FrameSize,1),1,SAllC,'same')),2);
    end
end 

Coinc =  sum(sign(conv2(ones(FrameSize,1),1,AllC,'same')),2);
p = zeros(size(Coinc));

% calculate the empirical p values times series for experimental and
% alternate coincidence time series.
AltP = zeros(size(AlternativeCoincs));
for i = 1:length(p)
    S = AlternativeCoincs(i,:);
    dist = sort(S);
    K = zeros(size(S,1),size(S,2)+1);
    R = [Coinc(i) S];
    for j = 0:max(R)
        if ~isempty(K(R==j))
            if isempty(find(dist < j+1,1))
                K(R==j) = 0;
            else
                K(R==j) = find(dist < j+1,1,'last')/Iter;
            end
        end
    end
    p(i) = K(1);
    AltP(i,:) = K(2:end);
    %dist = AlternativeCoincs(i,:);
    %p(i) = max([length(dist(dist>Coinc(i))),1])/Iter;
end

CoincSurprise = log10((1-p)./p);
CoincSurprise(CoincSurprise>3) = 3;
CoincSurprise(CoincSurprise<-3) = -3;

% Coinc = Coinc;
CoincRank = p;
pVal = empDist(Coinc,AlternativeCoincs);

NPCscore = -log10(pVal + (1/Iter));

% 	how to assess difference between the actual distribution of 
%   coincidental activity against that of the shifted alternatives?
% 	actually comparing the cummulative distributions of actlvl
% 		considered comparing cummulative distributions of actlvl ranks
% 		but it depends even more on the shuffle window. not to worry about for now.

altpVal = empDist(CoincRank,AltP); 
altNPCscore = -log10(altpVal + (1/Iter));
%%%%%

function p = empDist(empCoinc,empAlternatives)

% function p = empDist(Coinc,empLike)
% non-parametric statistic evaluating the distribution of an experimental
% distribution of coincidences (Coinc, a column time series of counts),
% against the distributions of alternative coincidence series (empLike), using rank
% euclidean distance from mean alternative coincidence distribution.
% designed to compliment outputs of localSurprise.m
% corrected to use cummulative distributions. As always intended

% Finn Upham 2016 02 18

m = max(max([empCoinc empAlternatives]));
x = 0:max(max([empCoinc empAlternatives]));
if m <= 1
    x = (0:25)/25;
end
rH = hist(empCoinc,x)';
cdistf = cumsum(rH)/sum(rH);
L = size(empAlternatives);
V = zeros(length(x),L(2));

for i = 1:L(2)
    c = hist(empAlternatives(:,i),x)'; %/L(1);
    V(:,i) = cumsum(c)/sum(c);
end

Vmean = mean(V,2);
dMA = sum((cdistf-Vmean).^2)^0.5;
dMV = sum((V-repmat(Vmean,1,L(2))).^2,1).^0.5;

p = length(dMV(dMV>dMA))/L(2);
return

