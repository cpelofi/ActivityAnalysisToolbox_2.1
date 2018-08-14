%% This script is a demo of the Activity Analysis toolbox
% first load the first example data set
% Finn Upham 2016 - 01 - 27
load data/Korhonen.mat

%% select one of the collections for analysis 
collN = 6;
D = Kor{collN};
P = D;

figure(1)
subplot(2,1,1)
% Plot of all ratings as Valence values (0.5 as crossing point for this
% bipolar scale) against time
plot(D.Time, D.Data)
title([D.Measure ' ratings of ' D.Piece ' in ' D.Audience])
xlabel('Time (s)'); ylabel('Rating range');
axis([D.Time(1) D.Time(end) 0 1])

subplot(2,1,2)
% Plot of all ratings as Valence values as colour scales in rows (heat map)
imagesc(D.Time,[1:D.Np],D.Data')
caxis([0 1])
title([D.Measure ' ratings of ' D.Piece ' in ' D.Audience])
xlabel('Time (s)'); ylabel('Response #');

figure(2)
subplot(3,1,1)
% Same collection of ratings, reporting the first order differences
plot(D.Time(2:end), diff(D.Data))
title([num2str(1/D.sF,3) ' Hz Differenced ' D.Measure ' ratings of ' D.Piece ' in ' D.Audience])
xlabel('Time (s)'); ylabel('Rating range');
diffRange = max(max(abs(diff(D.Data))));
axis([D.Time(1) D.Time(end) -diffRange diffRange])

subplot(3,1,2)
% each response in rows with valence scale (Colour) and time (X)
imagesc(D.Time(2,end),[1:D.Np],diff(D.Data)')
caxis([-0.1 0.1])
title([num2str(1/D.sF,3) ' Hz Differenced ' D.Measure ' ratings of ' D.Piece ' in ' D.Audience])
colorbar('East')
xlabel('Time (s)'); ylabel('Response #');

subplot(3,1,3)
% Activity-level times series for increases and decreases 
FrameSize = 2*D.sF;
HopSize = 2*D.sF;
Thresh = 0.025;
[AC,~,dT]=activityCount(D.Time,D.Data,FrameSize,HopSize,Thresh,'Inc');
AC(:,2)=activityCount(D.Time,D.Data,FrameSize,HopSize,Thresh,'Dec');
bar(dT,[-AC(:,2) AC(:,1)])
axis([D.Time(1) D.Time(end) -1 1])
set(gca,'yTick',[-1 -0.5 0 0.5 1],'yTickLabel',[1 0.5 0 0.5 1]) 
ylabel('Activity level');
xlabel('Time (s)')
title(['Increasing and Decreasing rating change activity levels in ' num2str(FrameSize/D.sF,3) ' s windows of synchrony'])
legend('Dec','Inc')
Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC(:,1),D.Np,Nbins);
C = coordScoreSimple(D.Time,D.Data,FrameSize,Thresh,'Inc',Nbins);
text(10,0.85,strcat('Inc C = ',num2str(C,3)));

[~,~,DAct]=simpleActivityTest(AC(:,2),D.Np,Nbins);
C = coordScoreSimple(D.Time,D.Data,FrameSize,Thresh,'Dec',Nbins);
text(10,-0.85,strcat('Dec C = ',num2str(C,3)));

C = coordScoreAlternating(D.Time,D.Data,FrameSize,Thresh,'Inc',Thresh,'Dec');
text(100,-0.85,strcat('Alt C = ',num2str(C,3)));

%%

figure

subplot(3,1,2)
FrameSize = 2*D.sF;
HopSize = 2*D.sF;
Thresh = 0.025;
[AC,~,dT]=activityCount(D.Time,D.Data,FrameSize,HopSize,Thresh,'Inc');
bar(dT,[AC(:,1)])
axis([D.Time(1) D.Time(end) 0 1])
ylabel('Activity level');
xlabel('Time (s)')
legend('Inc')
set(gca,'yTick',[0 0.5 1],'yTickLabel',[0 0.5 1]) 
title('Rating Increase of >2.5% in 2s time frames, non-overlapping')

subplot(3,1,3)
FrameSize = 2*D.sF;
HopSize = 2*D.sF;
Thresh = 0.025;
[acts,dT]=multiActPlotter(D.Time,D.Data,FrameSize,HopSize,'Inc',[0.1 0.05 0.025 0.01]);
area(dT,acts)
axis([D.Time(1) D.Time(end) 0 1])
ylabel('Activity level');
xlabel('Time (s)')
legend('Inc 0.1','Inc 0.05','Inc 0.025','Inc 0.01')
title('Rating Increase of many thresholds in 2s time frames, non-overlapping')

subplot(3,1,1)
plot(D.Time, D.Data)
hold on
plot(D.Time, mean(D.Data,2),'k.')
title([D.Measure ' ratings of ' D.Piece ' in ' D.Audience])
xlabel('Time (s)'); ylabel('Rating range');
axis([D.Time(1) D.Time(end) 0 1])

%% Demonstrating different kinds of activity evaluated by the activity count function

figure

subplot(3,1,1)
% ratings and their average time series
plot(D.Time, D.Data)
hold on
plot(D.Time, mean(D.Data,2),'k.')
title([D.Measure ' ratings of ' D.Piece ' in ' D.Audience])
xlabel('Time (s)'); ylabel('Rating range');
axis([D.Time(1) D.Time(end) 0 1])

subplot(3,1,3)
FrameSize = 2*D.sF;
HopSize = 2*D.sF;
Thresh = 0.025;
[acts,dT]=multiActPlotter(D.Time,D.Data,FrameSize,FrameSize,'Percent',[90 80 70 60 50]);
area(dT,acts);
axis([D.Time(1) D.Time(end) 0 1])
ylabel('Activity level');
xlabel('Time (s)')
legend('90th','80th','70th','60th','50th')
title('Rating values in the Nth percentil of used range, 2s time frames, non-overlapping')


subplot(3,1,2)
FrameSize = 2*D.sF;
HopSize = 2*D.sF;
Thresh = 0.025;
[acts,dT]=multiActPlotter(D.Time,D.Data,FrameSize,FrameSize,'UBound',[90 80 70 60 50]*0.01);
area(dT,acts);

axis([D.Time(1) D.Time(end) 0 1])
ylabel('Activity level');
xlabel('Time (s)')
legend('0.9','0.8','0.7','0.6','0.5')
title('Ratings above rating thresholds, 2s time frames, non-overlapping')




%% displaying all activity options
% note some of these are designed for signals other than ratings. 

options ={'Change','Inc','Dec','UBound','Xup','LBound','Percent'};
threshes = [0.025 0.025 0.025 0.5 0.5 0.5 90];
figure(3)
subplot(4,2,1)
respN = 5;
plot(D.Time, D.Data(:,respN))
axis([D.Time(1) D.Time(end) 0 1])
title('Single Response')

FrameSize = 2*D.sF;
HopSize = 2*D.sF;
Thresh = 0.025;
for i = 1:7
    subplot(4,2,1+i)
    [~,AllC,dT]=activityCount(D.Time,round(50*D.Data)/50,FrameSize,HopSize,threshes(i),options{i});
    stem(dT,AllC(:,respN));
    title(['Action ' options{i} ' Threshold ' num2str(threshes(i)) ' Framesize ' num2str(FrameSize/D.sF) 's'])
    axis([D.Time(1) D.Time(end) 0 1.1])
    set(gca,'yTick',[])
end

%% Activity tests and Coordination Score for simple activity of one type of event

P = Kor{1};
FrameSize = 2*P.sF;
HopSize = 2*P.sF;
Thresh = 0.025;

figure
spC = 3; spR = 1; spN = 1;

subplot(spC,spR,spN); spN = spN + 1;
[AC,~,dT]=activityCount(P.Time,P.Data,FrameSize,HopSize,Thresh,'Inc');
bar(dT,AC,0.5)
axis([dT(1) dT(end) 0 1])
ylabel('Activity level');
xlabel('Time (s)')
title([P.Measure ' rating increases activity-level times series in ' P.Piece ' excerpt'])

Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC,P.Np,Nbins);
subplot(spC,spR,spN); spN = spN + 1;
bar((0:P.Np)/P.Np,DAct,1)
axis([-0.02 1 0 0.35*length(AC)])
ylabel('# time frames')
xlabel('Activity Level')
title('Rating increases activity-level distributions, experimental and model')
legend Experimental Random
hold on
plot(((0:P.Np)-0.15)/P.Np,DAct(:,1),'.-')
plot(((0:P.Np)+0.15)/P.Np,DAct(:,2),'k.-')

spN = 2*(spC - 1)+1;
subplot(spC,2,spN); spN = spN + 1;
Nbins = 4;
[C,pval,DAct,bins]=simpleActivityTest(AC,P.Np,Nbins);
bar(1:Nbins,bins,1)
axis([0 Nbins+1 0 0.5*length(AC)])
hold on
plot([0 Nbins+1],[1 1]*length(AC)/Nbins,'k:')
ylabel('# time frames')
xlabel('GoF test bins')
title('GoF test bins 1')
%legend Experimental Random
text(Nbins-2,0.4*length(AC),['\chi^2 = ',num2str(C,3)]);
text(Nbins-2,0.35*length(AC),['d.f. = ' int2str(Nbins-2)]);
text(Nbins-2,0.3*length(AC),['p << 0.001 ']);
P.Time,P.Data,FrameSize,Thresh,'Inc'
C = coordScoreSimple(P.Time,P.Data,FrameSize,Thresh,'Inc',Nbins)
text(Nbins-2,0.45*length(AC),strcat('C Score (Inc) = ',num2str(C,3)));


subplot(spC,2,spN); spN = spN + 1;
Nbins = 5;
[C,pval,DAct,bins]=simpleActivityTest(AC,P.Np,Nbins);
bar(1:Nbins,bins,1)
hold on
plot([0 Nbins+1],[1 1]*length(AC)/Nbins,'k:')
axis([0 Nbins+1 0 0.5*length(AC)])
ylabel('# time frames')
xlabel('GoF test bins')
title('GoF test bins 2')
%legend Experimental Random
text(Nbins-2,0.4*length(AC),['\chi^2 = ' num2str(C,3)]);
text(Nbins-2,0.35*length(AC),['d.f. = ' int2str(Nbins-2)]);
text(Nbins-2,0.3*length(AC),['p << 0.001 ']);
C = coordScoreSimple(P.Time,P.Data,FrameSize,Thresh,'Inc',Nbins)
text(Nbins-2,0.45*length(AC),strcat('C Score (Inc) = ',num2str(C,3)));

%% Alternating activity test (not discussed in AA paper)

figure
spC = 5; spR = 1; spN = 1;

P = Kor{2};
FrameSize = 2*P.sF;
HopSize = 2*P.sF;
Thresh = 0.025;

FrameSize = 2; HopSize = 2; Thresh = 0.025;
[AC,AllC1,dT]=activityCount(P.Time,P.Data,FrameSize,HopSize,Thresh,'Inc');
[AC(:,2),AllC2,dT]=activityCount(P.Time,P.Data,FrameSize,HopSize,Thresh,'Dec');


% subplot 1: Rating Increases and Decreases activity-level times series
subplot(spC,spR,spN); spN = spN + 1;
bar(dT,[AC(:,1) -AC(:,2)])
axis([P.Time(1) P.Time(end) -1 1])
legend('Inc','Dec')
ylabel('Activity level');
xlabel('Time (s)')
title(['Increasing and Decreasing ' P.Measure  'rating change activity levels in one collection to ' P.Piece])

% subplot 2: Total rating change and Inc/Dec ratio time series. 
subplot(spC,spR,spN); spN = spN + 1;
bar(dT,[AC(:,1)+AC(:,2) -AC(:,1)./(AC(:,1)+AC(:,2))])
axis([P.Time(1) P.Time(end) -1 1])
legend('Total','AltRatio')
ylabel('Activity level');
xlabel('Time (s)')
title('Inc/Dec activity ratio and total rating change activity-levels from one collection')

% subplot Big Square to the left: this plots the concentration of time
% frames with specific Inc/Dec ratios and total rating change activity
% levels, from the actual collection's time series
spR = 5;
subplot(spC,spR,2*spR+[2:3 7:8])
[C,pval,DAct,Bins,v1,v2]=alternatingActivitiesTest(AllC1,AllC2),3;
g1 = [mean([v1{1}(end) v1{2}(1)]) mean([v1{2}(end) v1{3}(1)])];
g2 = [mean([v2{1}(end) v2{2}(1)]) mean([v2{2}(end) v2{3}(1)])];
[Co,h] = contourf(squeeze(DAct(:,:,2)),max(max(DAct(:,:,2))),':');
title('Actual Joint distribution')
a = caxis;
set(gca,'yTick',g1,'yTickLabel',[]) 
set(gca,'xTick',g2,'xTickLabel',[])
grid on

% subplot Big Square to the right: this plots the concentration of time
% frames with specific Inc/Dec ratios and total rating change activity
% levels, as expected if increases and decreases were not coordinated.

subplot(spC,spR,2*spR+[4:5 9:10])
[Co,h] = contourf(round(squeeze(DAct(:,:,1))),round(max(max(squeeze(DAct(:,:,1))))),':');
caxis(a)
colorbar
set(gca,'yTick',g1,'yTickLabel',[]) 
set(gca,'xTick',g2,'xTickLabel',[])
grid on
title('Random model joint distribution for total activity and ratio')
% colorbar
% set(h,'LineWidth',0);

% actual distribution of time freams by increases/decreases ratio
subplot(spC,spR,1+spR*[2:(spC-2)])
barh((0:size(AllC1,2))/size(AllC1,2),sum(DAct(:,:,1),2),1)
axis([0 0.2*length(AC) -0.02 1.02 ])
title('Distribution by Inc/Dec ratio')
ylabel('Inc/Dec activity level ratios')

subplot(spC,spR,spR*(spC-1)+(2:(spR-2)))
bar((0:size(AllC2,2))/size(AllC2,2),sum(DAct(:,:,1),1),1)
axis([-0.02 1.02 0 0.25*length(AC)])
title('Total activity levels (Inc and Dec) distribution')
xlabel('Total activity levels')
ylabel('# time frames')

subplot(spC,spR,spR*spC-1)
imagesc(Bins(:,:,1))
view([0 90])
set(gca,'yTick',[]) 
set(gca,'xTick',[]) 
axis xy
colorbar
title('Test table, Actual')
caxis([5 0.2*length(AC)])
%title('Contingency tables for Joint Distribution GoF test') 
subplot(spC,spR,spR*spC)
imagesc(Bins(:,:,2));
view([0 90])
set(gca,'yTick',[]) 
set(gca,'xTick',[]) 
axis xy
colorbar
title('Test table, Ind.')
caxis([5 0.2*length(AC)])

subplot(spC,spR,spR*(spC-1)+1)
CS = coordScoreAlternating(P.Time,P.Data,FrameSize,Thresh,'Inc',Thresh,'Dec',3);
text(0.2,0.9,['C score (Alt) = ' num2str(CS,3)]);
text(0.2,0.7,['\chi^2 = ' num2str(C,3)]);
text(0.2,0.5,['d.f. = ' int2str(4)]);
text(0.2,0.3,['pval = ' num2str(pval,3)]);
set(gca,'yTick',[]) 
set(gca,'xTick',[]) 

%% two Activities

% figure 5 same stim diff part Coll 1 ActLvl (inc) 1, Coll 2, ActLvl (Inc)
% 2) vert dist inc 1, middle square both, horiz inc 2. Column of together,
% alternating, and random examples.

figure
spC = 4; spR = 1; spN = 1;

P = Kor{1};
FrameSize = 2*P.sF;
HopSize = 2*P.sF;
Thresh = 0.025;

[AC,AllC1,dT]=activityCount(P.Time,P.Data,FrameSize,HopSize,Thresh,'Inc');
P = Kor{2};
[AC(:,2),AllC2,dT]=activityCount(P.Time,P.Data,FrameSize,HopSize,Thresh,'Inc');

subplot(spC,spR,spN); spN = spN + 1;
bar(dT,[AC(:,1) -AC(:,2)])
axis([P.Time(1) P.Time(end) -1 1])
ylabel('Activity level');
xlabel('Time (s)')
title('Increasing rating change activity levels for Arousal and Valence dimensions')
legend('Coll 1','Coll 2')
set(gca,'yTick',[-1 -0.5 0 0.5 1],'yTickLabel',[1 0.5 0 0.5 1]) 


spR = 4;
subplot(spC,spR,1+spR*[1:(spC-2)])
Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC(:,1),size(AllC1,2),Nbins);
barh((0:size(AllC1,2))/size(AllC1,2),DAct(:,1),1)
axis([0 0.4*length(AC) -0.02 1 ])
ylabel('Increases Activity levels')
xlabel('# frames')
title('Collection 1')

subplot(spC,spR,spR*(spC-1)+(2:(spR-1)))
Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC(:,2),size(AllC2,2),Nbins);
bar((0:size(AllC2,2))/size(AllC2,2),DAct(:,1),1)
axis([-0.02 1 0 0.4*length(AC)])
xlabel('Increases Activity levels')
ylabel('# frames')
title('Collection 2')

%subplot(spC,spR,spR+[2:4 7:9 12:14])
subplot(spC,spR,spR+[2:3 6:7])
Nbins = 3;
[C,p,DAct,Bins,v1,v2]=relatedActivitiesTest(AllC1,AllC2,Nbins);
g1 = [mean([v1{1}(end) v1{2}(1)]) mean([v1{2}(end) v1{3}(1)])];
g2 = [mean([v2{1}(end) v2{2}(1)]) mean([v2{2}(end) v2{3}(1)])];
[Co,h] = contourf(squeeze(DAct(:,:,2)),max(max(DAct(:,:,2))),':');
colorbar
title('Joint Activity distribution, Increases to Arousal and Valence dimensions')
set(gca,'yTick',g1,'yTickLabel',[]) 
set(gca,'xTick',g2,'xTickLabel',[])
grid on


%subplot(spC,spR,10)
subplot(spC,spR,spR*(spC-1)+1)
[Co,h] = contourf(squeeze(DAct(:,:,1)),round(max(max(DAct(:,:,1)))),':');
caxis([0 8])
axis xy
title('Independent model joint-Dist')
set(gca,'yTick',g1,'yTickLabel',[]) 
set(gca,'xTick',g2,'xTickLabel',[])
grid on
 
%subplot(spC,spR,15)
subplot(spC,spR,12)
imagesc(Bins(:,:,1))
caxis([0 20])
view([0 90])
set(gca,'yTick',[]) 
set(gca,'xTick',[]) 
axis xy
colorbar
title('Test table, Actual')

%subplot(spC,spR,20)
subplot(spC,spR,8)
imagesc(Bins(:,:,2));
caxis([0 20])
view([0 90])
axis tight
set(gca,'yTick',[]) 
set(gca,'xTick',[]) 
colorbar
title('Test table, Ind.')

subplot(spC,spR,spC*spR)
axis([0 1 0 1])
set(gca,'yTick',[]) 
set(gca,'xTick',[]) 
text(0.2,0.75,['\chi^2 = ' num2str(C,3)]);
text(0.2,0.5,['d.f. = ' int2str(4)]);
text(0.2,0.25,['pval = ' num2str(p,3)]);
P = Kor{1};
Q = Kor{2};
CS = coordScoreRelated(P.Time,P.Data,Q.Data,FrameSize,Thresh,'Inc',Thresh,'Inc',Nbins);
text(0.2,0.95,['C Score Ind = ' num2str(CS,3)]);
