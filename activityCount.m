function [AC,AllC,dT]=activityCount(Time,Series,FrameSize,HopSize,Thresh,option)

%function [AC,AllC,dT]=activityCount(Time,Series,FrameSize,HopSize,Thresh,option)
%
% function to identify the occurances and popularity of activity events 
% (specified by 'option' and 'Thresh') across the columns of "Series" over
% frames of size FrameSize + 1 sample at intervals of Hopsize. Columns are
% treated as signals to be translated into point processes.
%
% AC reports the proportion of columns showing the specified activity in
% each time frame.
% AllC reports the instances of activity events within each column/series 
% per timeframe
% dT reports the midpoint timestamp of each timeframe of AC and AllC.
%
% Action option must be one of:
% 'Change','Inc', 'Dec','UBound','LBound','Percent','LMax','LMin','Xup' 
% If the type of activity specified by option is found in a column to a
% degree exceeding the threshold 'thresh', that column is marked as active
% (1) for that time frame centered at the timestampr reported in the 
% corresponding row of the dT series. 
%
% 'Change' : whether absolute difference of series values from begining to end of
% timeframe is greater than Thresh
% 'Inc' : whether difference of series values from begining to end of
% timeframe is greater than Thresh
% 'Dec' : whether difference of series values from begining to end of
% timeframe is less than -Thresh
% 'UBound' : whether series values exceeds Thresh for at least one sample
% in timeframe
% 'LBound' : whether series values falls below Thresh for at least one sample
% in timeframe
% 'Percent': whether average value in time frame exceeds (or falls below
% for Thresh < 0 ) the percentile of each column's distribution of values
% 'LMax' : whether there is a local maxima with signal values above thresh
% 'LMin' : whether there is a local minima with signal values above thresh
% 'Xup' : whether the signal rises across threshole in the timeframe

% Version 4.0 (earlier versions are compatible with older scripts)
% Finn Upham 2016 01 20
 
if nargin==5
    option='Change'; % change is the default activity type
end

if size(Time,1)~=size(Series,1)
error('Matlab:actionCount:RowsMatch',...
        'Rows of Time and Series must match')
end

% set up the frames for activity analysis and define the time points
% associated with each frame.

k=round(HopSize);
l=round(FrameSize);

if k>l
error('Matlab:actionCount2:HopSize',...
        'Warning hopsize larger than framesize, may miss activity.')
end

T=Time(round(k/2):k:size(Series,1),:);
dT=T;
frame=cell(length(T),1);
cSize = size(Series,2);

for i=1:length(T)
    if k*(i-1)+l>size(Series,1)
        frame{i}=Series(k*(i-1):end,:);
    elseif i ==1
        frame{i}=Series(1:l,:);
    else
        frame{i}=Series(k*(i-1):k*(i-1)+l,:);
    end
end

%Now evaluate and save activity over each frame, according to the requested
%activity type

AC=zeros(size(T));
AllC=zeros(size(T,1),size(Series,2));

switch option
    case 'Change'
        % whether absolute difference of series values from begining to end 
        % of timeframe is greater than Thresh
        for i=1:length(T)
            %take the difference of the first and last frame    
            D=abs(frame{i}(end,:)-frame{i}(1,:));
            %Set D(r)==1 if response r changes over frame i, else D(r)==0
            D(D<Thresh)=0;
            D=abs(sign(D));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'Inc'
        % 'Inc' : whether difference of series values from begining to end 
        % of timeframe is greater than Thresh
        for i=1:length(T)
            %take the difference of the first and last frame    
            D=frame{i}(end,:)-frame{i}(1,:);
            % Set D(r)==1 if response r increased over frame i at least Thresh,
            % else D(r)==0
            D(D<Thresh)=0;
            D=sign(D);
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'Dec'
        % 'Dec' : whether difference of series values from begining to end 
        % of timeframe is less than -Thresh
        for i=1:length(T)
            %take the difference of the first and last frame    
            D=frame{i}(end,:)-frame{i}(1,:);
            %Set D(r)==1 if response r decreased over frame i, else D(r)==0
            D(D>-Thresh)=0;
            D=abs(sign(D));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'UBound'
        % 'UBound' : whether series values exceeds Thresh for at least one
        % sample in timeframe
        for i=1:length(T)
            %Set the UBound to zero
            D=frame{i};
            D=D-Thresh;
            %Push all non-negative response values to 1, and rest to 0
            D(D>=0)=1;
            D(D<0)=0;
            %sum over frame to catch which response have any points with values
            %over or at the UBound and set their activity values to 1 for the frame
            D=sign(sum(D,1));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'LBound'
        % 'LBound' : whether series values falls below Thresh for at least
        % one sample in timeframe
        for i=1:length(T)
            %Set the LBound to zero
            D=frame{i};
            D=D-Thresh;
            %Push all non-positive response values to 1, and rest to 0
            D(D<=0)=0;
            D(D>0)=1;
            D = 1-D;
            %sum over frame to catch which response have any points with values
            %below or at the LBound and set their activity values to 1 for the frame
            D=sign(sum(D,1));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'Percent'
        %to identify the frames with the highest average values (in place of
        %sample wise rank value), we first average response values over each
        %frame
        V=AllC;
        for i=1:length(T)
            V(i,:)=mean(frame{i});
        end
        if Thresh<0
            V=-V;
            Thresh=-Thresh;
        end
        %find the percentil values for each response over these frame averages
        Y=prctile(V,Thresh);
        %throwing out instances when the percentile includes too much of the
        %signal
        for i=1:size(V,2)
            if length(V(V(:,i)>=Y(i)))>length(T)*2*(100-Thresh)/100
                V(:,i)=Y(i)-1;
            end
        end
        for i=1:length(T)
            %per frame, treat the percentiles values as UBounds per responses
            D=V(i,:)-Y;
            %Push all non-negative response values to 1, and rest to 0
            D(D>=0)=1;
            D(D<0)=0;
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'LMax' 
        % pick local maxima (peaks) which crest above thresh in series units.
        for i=1:length(T)
            %Set the LBound to zero
            if k<l
                D=frame{i};
            elseif k==l
                if i < length(T)
                    D=[frame{i};frame{i+1}(1:2,:)];
                end
            end
            D2 = D;
            D = diff([D; D(end,:)],1,1);
            D = diff(sign([D(1,:);D]),1,1);
            D(end,:) = 0;
            D(D>0) = 0;
            D(D2<Thresh) = 0;
            D = abs(sign(sum(D,1)));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'LMin'  
        % pick local minima (valleys) which dip below thresh in series units.
        for i=1:length(T)
            %Set the LBound to zero
            if k<l
                D=frame{i};
            elseif k==l
                if i < length(T)
                    D=[frame{i};frame{i+1}(1:2,:)];
                end
            end
            D2 = D;
            D = diff([D; D(end,:)],1,1);
            D = diff(sign([D(1,:);D]),1,1);
            D(end,:) = 0;
            D(D<0)=0;
            D(D2>Thresh)=0;
            D=abs(sign(sum(D,1)));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    case 'Xup'
        % whether the signal rises across threshold in the timeframe
        for i=1:length(T)
            if k<l
                D=frame{i};
            elseif k==l
                if i < length(T)
                    D=[frame{i};frame{i+1}(1,:)];
                end
            end
            %Set the threshold to zero to capture crossing upwards
            D = sign(D-Thresh);
            D=diff(D,1,1);
            D(D<0)=0;
            D=abs(sign(sum(D,1)));
            %Save the individual activity and activity count for frame i
            AllC(i,:)=D;
            AC(i)=sum(D)/cSize;
        end
    otherwise
         error('MATLAB:actionCount:ppOutput', ...
                'option not supported')
end




