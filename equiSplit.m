function [v,d] = equiSplit(V,n,minim)

% function [v,d] = equiSplit(V,n,minim)

% V is a histogram output, an integer-valued list, 
% n is the target number of bins the histogram should be simplified to
% minim is minimum number of samples per aggregate bin  
% (e.g. 5 for Goodness of Fit test).
% 
% This function splits V into n intervals with maximally even values per
% interval, all exceeding min. If a division is not possible given the min
% provided, a warning reported and the histogram is instead sorted into 
% (n-1) bins. 
% v is a cell with size equal to the number of bins, each entry the indecies 
% of V for a single bin. If no solution satisfied the minimum
% bin size using 2 bins, v is returned empty.

% Finn Upham, August 22nd, 2012

% make V a column vector
if size(V,1) == 1
    V = V';
end
if size(V,2) > 1
    error('V must be a column vector')
end

% generate the cummulative distribution of V
cumV = cumsum(V);
N = sum(V);
L = length(V);

% Select 2 indices close to the even distribution of values for each break
% point.
A = repmat(cumV,1,n-1) - repmat(N*(1:n-1)/n,L,1);
[~,Bind] = sort(abs(A),1,'ascend');
Bind = Bind(1:2,:);

% Find combination of break points which 1) satisfy minim 2)spread most
% evenly across the bins.
binOption = zeros(n+1,2^(n-1));
ind = zeros(n-1,1);

for i = 0:2^(n-1)-1
    b = ind2binV(i,n-1)+1;
    for j = 1:n-1
        ind(j) = Bind(b(j),j);
    end
    binOption(:,i+1) = [0;cumV(ind);N];
end

options = 1:2^(n-1);
% check which satisfy min

minCheck = min(diff(binOption),[],1);
options = options(minCheck>=minim);

% if no options survive this criteria, repeat process with smaller n
if isempty(options)
    if n > 2
        [v,d] = equiSplit(V,n-1,minim);
    else
        v = [];
        d = cell(0);
        %fprintf('\nNo cutting this distribution.\n');
    end
end


% if there are options, choose the most even set
% hun. I've been using N1 for this. Maybe should use N2.
if ~isempty(options)
    %diffCheck = sum(abs(diff(binOption(:,options))));
    diffCheck = sum(diff(binOption(:,options)).^2);
    [~,minDiff] = min(diffCheck);
    bestOption = options(minDiff);
    b = ind2binV(bestOption-1,n-1)+1;
    for j = 1:n-1
        ind(j) = Bind(b(j),j);
    end
    ind = [0;ind;L];
    v = cell(n,1);
    d = zeros(n,1);
    for i = 1:n
        v{i} = ind(i)+1:ind(i+1);
        d(i) = sum(V(v{i}));
    end

end

return

function b = ind2binV(ind,N)
% ind2binV converts a single integer (index) into a vector of 
    binStr = dec2bin(ind,N);
    
    b = zeros(1,N);
    
    for bin = 1:length(binStr)
        b(bin) = str2double(binStr(bin));
    end
return
