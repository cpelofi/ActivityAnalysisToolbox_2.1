% this script requires a data set which can be downloaded from this link:
% http://www.mediafire.com/download/tuzwtaa94bo7a0z

fprintf('This script uses a data set which must be downloaded seperately.\n http://www.mediafire.com/download/tuzwtaa94bo7a0z\n')

clear
load data/stim4.mat

fprintf('This data reports the continuous felt emotion ratings \n and some psychophysiological measures recorded during\n repeated listenings by one participant to\n The Be Good Tanya''s song The Littlest Birds \n Youtube link: https://youtu.be/VdIhpkEkC4c\n')

fprintf('Analysis of psychophysiological signals follow criteria\n articulated in ICMPC presentation:\n   Demo: Activity Analysis on Psychophysiological\n   Measures of Response to Music.\n Video and slides: https://wp.me/ptLMM-PL\n')

%% Basic summary plots

time = stim.time;
track = stim.track;
kStep = stim.time.Hz/5; % simple downsample to 5 Hz

beh = stim.beh.aro;

figure(1)
clf
subplot(4,1,1)
plot(time.s,beh.D)
axis([time.s(1) time.s(end) 0 1])
title([track.composer ' ' track.piece ' ' beh.T ' ratings and average'])
hold on 
plot(time.s,mean(beh.D,2),'k.')
xlabel(time.sT)
ylabel(beh.T)

subplot(4,1,2)
imagesc(time.s(1:kStep:end,:),[1:size(beh.D,2)],beh.D(1:kStep:end,:)')
caxis([0 1])
colorbar
title([track.composer ' ' track.piece ' ' beh.T ' ratings, each session'])
xlabel(time.sT)
ylabel('Session')

beh = stim.beh.val;
subplot(4,1,3)
plot(time.s,beh.D)
axis([time.s(1) time.s(end) 0 1])
title([track.composer ' ' track.piece ' ' beh.T ' ratings and average'])
hold on 
plot(time.s,mean(beh.D,2),'k.')
xlabel(time.sT)
ylabel(beh.T)

subplot(4,1,4)
imagesc(time.s(1:kStep:end,:),[1:size(beh.D,2)],beh.D(1:kStep:end,:)')
caxis([0 1])
colorbar
title([track.composer ' ' track.piece ' ' beh.T ' ratings, each session'])
xlabel(time.sT)
ylabel('Session')

%% rating change coordination


time = stim.time;
track = stim.track;
kStep = stim.time.Hz/5; % simple downsample to 5 Hz
beh = stim.beh.aro;
beh.D = beh.D(1:kStep:end,:);
time.s = time.s(1:kStep:end,:);
time.Hz = 5;

figure

spC = 3; spR = 1; spN = 1;
subplot(spC,spR,spN); spN = spN + 1;
plot(time.s,beh.D)
hold on;
plot(time.s,mean(beh.D,2),'k','linewidth',2)
axis([time.s(1) time.s(end) 0 1])
ylabel('Rating Range')
xlabel('Time (s)')
title([track.composer ' ' track.piece ' ' beh.T ' ratings'])

% C 2 activity level time series, with results of joint Alt
subplot(spC,spR,spN); spN = spN + 1;
FrameSize = 2*time.Hz; HopSize = 2*time.Hz; Thresh = 0.025;
[AC,AllC1,dT]=activityCount(time.s,beh.D,FrameSize,HopSize,Thresh,'Inc');
[AC(:,2),AllC2]=activityCount(time.s,beh.D,FrameSize,HopSize,Thresh,'Dec');
bar(dT,[AC(:,1) -AC(:,2)])
axis([time.s(1) time.s(end) -1 1])
ylabel('Activity level');
xlabel('Time (s)')
title('Increasing and Decreasing rating change activity level times series')

% C 3 mono Inc, mono Dist, Joint dist
spR = 3; spN = (spC-1)*spR +1;
Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC(:,1),size(beh.D,2),Nbins);
subplot(spC,spR,spN); spN = spN + 1;
bar((0:size(beh.D,2))/size(beh.D,2),DAct/length(AC),1,'EdgeColor','None')
axis([-0.02 1 0 0.35])
ylabel('Frequency')
xlabel('Activity Level')
title('Increases activity level distributions')
legend Actual Random
hold on
plot(((0:size(beh.D,2))-0.15)/size(beh.D,2),DAct(:,1)/length(AC),'.-')
plot(((0:size(beh.D,2))+0.15)/size(beh.D,2),DAct(:,2)/length(AC),'k.-')
C = coordScoreSimple(time.s,beh.D,FrameSize,Thresh,'Inc',4);
text(0.4,0.15,strcat('C score = ',num2str(C,3)));

[~,~,DAct]=simpleActivityTest(AC(:,2),size(beh.D,2),Nbins);
subplot(spC,spR,spN); spN = spN + 1;
bar((0:size(beh.D,2))/size(beh.D,2),DAct/length(AC),1,'EdgeColor','None')
axis([-0.02 1 0 0.35])
ylabel('Frequency')
xlabel('Activity Level')
title('Decreases activity level distributions')
legend Actual Random
hold on
plot(((0:size(beh.D,2))-0.15)/size(beh.D,2),DAct(:,1)/length(AC),'.-')
plot(((0:size(beh.D,2))+0.15)/size(beh.D,2),DAct(:,2)/length(AC),'k.-')
C = coordScoreSimple(time.s,beh.D,FrameSize,Thresh,'Dec',4);
text(0.4,0.15,strcat('C score = ',num2str(C,3)));

subplot(spC,spR,spN); spN = spN + 1;
[~,~,DAct] = alternatingActivitiesTest(AllC1,AllC2,3);
[Co,h] = contourf(squeeze(DAct(:,:,2)),max(max(DAct(:,:,2))),':');
title('Actual Joint distribution')
caxis([0 max(max(DAct(:,:,2)))])
ylabel('Inc/Dec ratio')
xlabel('Total Activity Level')
set(gca,'yTick',[])
set(gca,'xTick',[])
C = coordScoreAlternating(time.s,beh.D,FrameSize,Thresh,'Inc',Thresh,'Dec');
text(14,14,strcat('C score = ',num2str(C,3)));

%% return to figure 1 to add activity analysis results on rating changes

time = stim.time;
track = stim.track;
kStep = stim.time.Hz/5; % simple downsample to 5 Hz
beh = stim.beh.aro;
beh.D = beh.D(1:kStep:end,:);
time.s = time.s(1:kStep:end,:);
time.Hz = 5;


figure(1)
clf
subplot(6,1,1)
plot(time.s,beh.D)
axis([time.s(1) time.s(end) 0 1])
title([track.composer ' ' track.piece ' ' beh.T ' ratings and average'])
hold on 
plot(time.s,mean(beh.D,2),'k.')
%xlabel(time.sT)
ylabel(beh.T)

subplot(6,1,2)
imagesc(time.s,[1:size(beh.D,2)],beh.D')
caxis([0 1])
colorbar
title([track.composer ' ' track.piece ' ' beh.T ' ratings, each session'])
%xlabel(time.sT)
ylabel('Session')

subplot(6,1,3)
FrameSize = 2*time.Hz; HopSize = 2*time.Hz; Thresh = 0.025;
[AC,AllC1,dT]=activityCount(time.s,beh.D,FrameSize,HopSize,Thresh,'Inc');
[AC(:,2),AllC2]=activityCount(time.s,beh.D,FrameSize,HopSize,Thresh,'Dec');
bar(dT,[AC(:,1) -AC(:,2)])
axis([time.s(1) time.s(end) -1 1])
ylabel('Activity level');
%xlabel('Time (s)')
title('Increasing and Decreasing rating change activity level times series')

Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC(:,1),size(beh.D,2),Nbins);
C = coordScoreSimple(time.s,beh.D,FrameSize,Thresh,'Inc',Nbins);
text(10,0.5,strcat('Inc C = ',num2str(C,3)));

[~,~,DAct]=simpleActivityTest(AC(:,2),size(beh.D,2),Nbins);
C = coordScoreSimple(time.s,beh.D,FrameSize,Thresh,'Dec',Nbins);
text(10,-0.5,strcat('Dec C = ',num2str(C,3)));

C = coordScoreAlternating(time.s,beh.D,FrameSize,Thresh,'Inc',Thresh,'Dec',3);
text(100,-0.75,strcat('Alt C = ',num2str(C,3)));



beh = stim.beh.val;
beh.D = beh.D(1:kStep:end,:);
time.Hz = 5;

subplot(6,1,4)
plot(time.s,beh.D)
axis([time.s(1) time.s(end) 0 1])
title([track.composer ' ' track.piece ' ' beh.T ' ratings and average'])
hold on 
plot(time.s,mean(beh.D,2),'k.')
%xlabel(time.sT)
ylabel(beh.T)

subplot(6,1,5)
imagesc(time.s,[1:size(beh.D,2)],beh.D')
caxis([0 1])
colorbar
title([track.composer ' ' track.piece ' ' beh.T ' ratings, each session'])
%xlabel(time.sT)
ylabel('Session')

subplot(6,1,6)
FrameSize = 2*time.Hz; HopSize = 2*time.Hz; Thresh = 0.025;
[AC,AllC1,dT]=activityCount(time.s,beh.D,FrameSize,HopSize,Thresh,'Inc');
[AC(:,2),AllC2]=activityCount(time.s,beh.D,FrameSize,HopSize,Thresh,'Dec');
bar(dT,[AC(:,1) -AC(:,2)])
axis([time.s(1) time.s(end) -1 1])
ylabel('Activity level');
%xlabel('Time (s)')
title('Increasing and Decreasing rating change activity level times series')
Nbins = 4;
[~,~,DAct]=simpleActivityTest(AC(:,1),size(beh.D,2),Nbins);
C = coordScoreSimple(time.s,beh.D,FrameSize,Thresh,'Inc',4);
text(10,0.5,strcat('Inc C = ',num2str(C,3)));

[~,~,DAct]=simpleActivityTest(AC(:,2),size(beh.D,2),Nbins);
C = coordScoreSimple(time.s,beh.D,FrameSize,Thresh,'Dec',4);
text(10,-0.5,strcat('Dec C = ',num2str(C,3)));

C = coordScoreAlternating(time.s,beh.D,FrameSize,Thresh,'Inc',Thresh,'Dec',3);
text(100,-0.75,strcat('Alt C = ',num2str(C,3)));

%% Basic summary plots

a = fieldnames(stim);
time = stim.time;
kStep = 5;

for i = 1:7
    eval(['b = fieldnames(stim.' a{i+3} ');'])
    figure(i+3)
    for j = 1:length(b)
        eval(['bio = stim. ' a{i+3} '.' b{j} ';'])
        subplot(length(b),2,2*(j-1)+1)
        plot(time.s,bio.D)
        title(['Raw Physio ' bio.T])
        axis tight
    
        subplot(length(b),2,2*j)
        imagesc(time.s(1:kStep:end,:),[1:size(bio.D,2)],bio.D(1:kStep:end,:)')
        title(['Raw Physio ' bio.T])
        colorbar
    end
end

%% Expressive emg


Np = 24;
kStep = 10;

figure
subplot(3,1,1)
bio = stim.emgc.rms;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Raw Physio ' bio.T])
ylabel('Session')
caxis([0 4])
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(3,1,2)
bio = stim.emgz.rms;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Raw Physio ' bio.T])
ylabel('Session')
caxis([0 4])
colorbar
hold on
plot([time.s(1) time.s(end)],[5 5],'k')
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(3,1,3)
bio = stim.emgt.rms;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Raw Physio ' bio.T])
ylabel('Session')
xlabel('Time (s)')
caxis([0 8])
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

%% local surprise

time = stim.time;
kStep = 5;
sF = 20;
Np = 24;
kStep = 10;

figure
subplot(5,1,1)
bio = stim.resp.filt;
data = diff(bio.D(1:kStep:end,:))';
imagesc(time.s(kStep:kStep:end,:),[1:Np],data)
title(['Differenced Physio ' bio.T ' over first 50s of stimulus'])
range = [prctile(prctile(data',1),5) prctile(prctile(data',99),75)];
caxis(range)
%colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')
axis([-10 50 0 24])

subplot(5,1,2)
bio = stim.resp.Inha;
[~,D,T] = activityCount(time.s,bio.D,100/sF,100/sF,0.01,'UBound');
imagesc(T,[1:Np],D'); caxis([0 2])
axis([-10 50 0 24])
title('Inspiration onsets')

subplot(5,1,3)
Time = T;
alpha = 0.025;
thresh = log10((1-alpha)/alpha);
sWindow = 0.3*sF;
bWindow = 10*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Inspiration onset activity levels (0.3 s,10 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(NPCscore)])
axis([-10 50 0 1])

subplot(5,1,4)
plot(Time,pCoinc,'b','linewidth',1)
axis([-10 50 0 1])
title(['Local pvalues of inspiration onset coincident (0.3 s,10 s), p = ' num2str(p)])

subplot(5,1,5)
plot(Time,sCoinc,'r','linewidth',1)
axis([-10 50 -3 3])
hold on
plot([-10 time.s(end)],[2 2],':k')
plot([-10 time.s(end)],[-2 -2],':k')
title(['Local Surprise of inspiration onset coincident (0.3 s,10 s), p = ' num2str(p)])
xlabel('Time(s)')

%% Synchronous activity trap and resp
%load stim5

a = fieldnames(stim);
time = stim.time;
Np = 24;
kStep = 10;
figure
subplot(4,1,1)
bio = stim.resp.filt;
data = diff(bio.D(1:kStep:end,:))';
imagesc(time.s(kStep:kStep:end,:),[1:Np],data)
title(['Differenced Physio ' bio.T])
range = [prctile(prctile(data',1),25) prctile(prctile(data',99),75)];
caxis(range)
%colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,2)
bio = stim.resp.Inha;
kStep = 5;
sF = 20;
[~,D,T] = activityCount(time.s,bio.D,100/sF,100/sF,0.01,'UBound');
Time = T;
sWindow = 0.3*sF;
bWindow = 10*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Inspiration onset activity levels (0.3 s,10 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(NPCscore,2)])
axis([-10 time.s(end) 0 1])

subplot(4,1,3)
bio = stim.emgt.filtN;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Physio ' bio.T])
caxis([0 2])
%colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,4)
[~,D,T] = activityCount(time.s,bio.D,5,5,2,'UBound');
Time = T;
sWindow = 1*sF;
bWindow = 20*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Trap contraction peak activity levels (1 s,20 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(NPCscore,2)])
axis([-10 time.s(end) 0 0.75])
xlabel('Time (s)')

%% phys else 

a = fieldnames(stim);
time = stim.time;
Np = 24;
kStep = 20;
figure
subplot(4,1,1)
bio = stim.bvp.bpm;
imagesc(time.s(kStep:kStep:end,:),[1:Np],diff(bio.D(1:kStep:end,:))')
title(['Differenced Physio ' bio.T])
caxis([-2 2])
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,2)
bio = stim.scon.filt;
data = diff(bio.D(1:kStep:end,:))';
imagesc(time.s(kStep:kStep:end,:),[1:Np],data)
title(['Differenced Physio ' bio.T])
range = [prctile(prctile(data',20),25) prctile(prctile(data',90),75)];
caxis(range)
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,3)
bio = stim.temp.filt;
imagesc(time.s(kStep:kStep:end,:),[1:Np],diff(bio.D(1:kStep:end,:))')
title(['Differenced Physio ' bio.T])
data = diff(bio.D(1:kStep:end,:))';
range = [prctile(prctile(data',1),25) prctile(prctile(data',99),75)];
caxis(range)
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,4)
bio = stim.resp.bpm;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Physio ' bio.T])
caxis([10 30])
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')
xlabel('Time (s)')


%% Synchronous activity Facial sEMG, Corrugator and Zygomaticus
%load stim5

% Respiration rate
a = fieldnames(stim);
time = stim.time;
Np = 24;
kStep = 20;
figure
subplot(4,1,1)
bio = stim.resp.bpm;
imagesc(time.s(kStep:kStep:end,:),[1:Np],diff(bio.D(1:kStep:end,:))')
title(['Differenced Physio ' bio.T])
caxis([-1 1])
%colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,2)
kStep = 20;
sF = 5;
[~,D,T] = activityCount(time.s(1:kStep:end,:),bio.D(1:kStep:end,:),2*sF,1,0.5,'Inc');
Time = T;
sWindow = 0.4*sF;
bWindow = 50*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Respiration rate increases activity levels (0.4 s,50 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(NPCscore,2)])
axis([-10 time.s(end) 0 1])

% Heart rate
kStep = 20;
subplot(4,1,3)
bio = stim.bvp.bpm;
imagesc(time.s(kStep:kStep:end,:),[1:Np],diff(bio.D(1:kStep:end,:))')
title(['Differenced Physio ' bio.T])
caxis([-5 5])
%colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,4)
kStep = 20;
sF = 5;
[~,D,T] = activityCount(time.s(1:kStep:end,:),bio.D(1:kStep:end,:),2*sF,1,4,'Dec');
Time = T;
sWindow = 2*sF;
bWindow = 40*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Heart rate decreases activity levels (2 s,40 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(NPCscore,2)])
axis([-10 time.s(end) 0 1])

%% Synchronous activity Skin Conductance and Finger Temperature
%load stim5

% Skin Conductance
a = fieldnames(stim);
time = stim.time;
Np = 24;
kStep = 10;
figure
subplot(4,1,1)
bio = stim.scon.filt;
data = diff(bio.D(1:kStep:end,:))';
imagesc(time.s(kStep:kStep:end,:),[1:Np],data)
title(['Differenced Physio ' bio.T])
range = [prctile(prctile(data',20),25) prctile(prctile(data',90),75)];
caxis(range)
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,2)
kStep = 20;
sF = 5;
bio = stim.scon.ORt;
[~,D,T] = activityCount(time.s,bio.D,2*100,100/sF,0.05,'UBound');
Time = T;
sWindow = 2*sF;
bWindow = 30*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(Time,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Orienting Responses activity levels (2 s,30 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(altNPCscore,2)])
axis([-10 time.s(end) 0 1])
%
% Heart rate
kStep = 20;
subplot(4,1,3)
bio = stim.temp.filt;
imagesc(time.s(kStep:kStep:end,:),[1:Np],diff(bio.D(1:kStep:end,:))')
title(['Differenced Physio ' bio.T])
data = diff(bio.D(1:kStep:end,:))';
range = [prctile(prctile(data',1),25) prctile(prctile(data',99),75)];
caxis(range)
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')


subplot(4,1,4)
kStep = 20;
sF = 5;
[~,D,T] = activityCount(time.s(1:kStep:end,:),bio.D(1:kStep:end,:),2*sF,1,0.01,'Dec');
Time = T;
sWindow = 1*sF;
bWindow = 40*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Finger Teperature decreases activity levels (2 s,40 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(altNPCscore,2)])
axis([-10 time.s(end) 0 1])


%% Synchronous activity Facial EMG, Corrugator, Zygomaticus
%load stim5

% Respiration rate
a = fieldnames(stim);
time = stim.time;
Np = 24;
kStep = 10;
figure
subplot(4,1,1)
bio = stim.emgc.rms;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Raw Physio ' bio.T])
ylabel('Session')
caxis([0 4])
colorbar
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,2)
kStep = 20;
sF = 10;
[~,D,T] = activityCount(time.s(1:kStep:end,:),bio.D(1:kStep:end,:),1,1,3.5,'UBound');
Time = T;
sWindow = 0.5*sF;
bWindow = 40*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Corrugator peaks activity levels (0.5 s, 40 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(altNPCscore,2)])
axis([-10 time.s(end) 0 1])

% Zygomaticus
kStep = 10;
subplot(4,1,3)
bio = stim.emgz.rms;
imagesc(time.s(1:kStep:end,:),[1:Np],bio.D(1:kStep:end,:)')
title(['Raw Physio ' bio.T])
ylabel('Session')
caxis([0 4])
colorbar
hold on
plot([time.s(1) time.s(end)],[5 5],'k')
hold on
plot([time.s(1) time.s(end)],[5:5:20;5:5:20]'+0.5,'w')

subplot(4,1,4)
kStep = 20;
sF = 10;
[~,D,T] = activityCount(time.s(1:kStep:end,:),bio.D(1:kStep:end,:),1,1,3.5,'UBound');
Time = T;
sWindow = 1*sF;
bWindow = 40*sF;
[p,Coinc,pCoinc,sCoinc,AltC,AltP,altpVal,NPCscore,altNPCscore]=localActivityTest(D,sWindow,bWindow,1000);
hold on
p975 = prctile(AltC',97.5)';
plot(T,p975/size(D,2),'b')
bar(Time,Coinc/size(D,2))
title(['Zygomaticus peaks activity levels (1 s, 40 s) and expected 97.5 percentile rank, p = ' num2str(p) ', C score = ' num2str(altNPCscore,2)])
axis([-10 time.s(end) 0 1])

