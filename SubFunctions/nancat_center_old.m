function [centData,cIdx,cVal,R,dR,maxR] = nancat_center(data,center,dim,Cent)
%% nancat_center: centers input data around given value and pads with Nan's so all vectors are the same size
%   INPUTS:
%       data       	:   input data (matrix or cell of vectors)
%       center     	:   data center value (set to empty when using Cent)
%       dim         :   dimension to center around
%       Cent     	:   center index , overrides center value
%   OUTPUTS:
%       centDara    :   centered data
%       cIdx        :   center indicies
%       cVal        :   computed center values used
%       R           :   indicies from center to ends of each column
%       dR         	:   Nan indicies added above and below each column
%---------------------------------------------------------------------------------------------------------------------------------
if nargin<3
    dim = 1; % defaults to columns
    if nargin<2
        center = 0; % default center value is 0
    end
elseif nargin==4
    if ~isempty(center)
        error('Error: the input "center" must be empty when using the input "Cent"')
    end
    center = 0;
end

if iscell(data) % for cell array
    if dim~=1 && dim~=2 % only works on first two dimesions
        error('Error: dimension must be 1 or 2')
    else
        if dim==2 % transform to column vectors if row vectors are specified
            data = data';
        end
        sizeCheck = cellfun(@(x) isvector(x), data, 'UniformOutput', true);
        
        if any(~sizeCheck,'all')
            Error('Error: arrays in cells must be vectors ')
        end
        
        data = cellfun(@(x) x(:), data, 'UniformOutput', false); % transform to column vectors if row vectors are specified
    end
    [n_data, n_array] = size(data);
    [n_center, ~] = cellfun(@(x) size(x), data, 'UniformOutput', true); % size of dimension to center around & other dimension
    
    if n_center<=2
        error('Error: need at least 3 values in array')
    end

    [cVal, ~] = cellfun(@(x) min(abs(x - center)), data, 'UniformOutput', true); % get value closest to center value
    R = cell(1,n_array); % store the # of indicies above and below the center value
    cIdx = nan(size(data));
    for jj = 1:n_array % for each column
        R{jj} = nan(2,n_data);
        for kk = 1:size(n_center,1)
            [cIdx(kk,jj),~] = find( data{kk,jj}==cVal(kk,jj) ); % find location of center value
            if nargin==4
                if isnan(Cent)
                    cIdx(kk,jj) = n_center(kk,jj);
                else
                    cIdx(kk,jj) = Cent;
                end
            end
            R{jj}(1,kk) = length(1:cIdx(kk,jj)-1); % length of data below center
            R{jj}(2,kk) = length((cIdx(kk,jj)+1):n_center(kk,jj)); % length of data above center
        end
    end

    dR = cellfun(@(x) diff(x), R, 'UniformOutput', false); % change in length needed to center data for each column
    maxR = cellfun(@(x) max(x,[],2), R, 'UniformOutput', false);
    
    n_new = cellfun(@(x) sum(x)+1, maxR, 'UniformOutput', true); % new length
    centData = cell(1,n_array); % centered data
    for jj = 1:n_array % for each column
    	centData{jj} = nan(n_new(jj),n_data); % centered data
        for kk = 1:n_data
            dR{jj,kk} = abs(maxR{jj} - R{jj}(:,kk));
            centData{jj}(:,kk) = cat_pad(data{kk,jj},dR{jj,kk},nan);
        end
    end
    
    if length(centData)==1
        centData = centData{1};
    end
    
elseif ismatrix(data) % for matrix
    if dim~=1 && dim~=2 % only works on first two dimesions
        error('Error: dimension must be 1 or 2')
    else
        if dim==2 % transform to column vectors if row vectors are specified
            data = data';
        end
    end
    [n_center, n_array] = size(data); % size of dimension to center around & other dimension

    if n_center<=2
        error('Error: need at least 3 values in array')
    end

    [cVal, ~] = min(abs(data - center)); % get value closest to center value
    R = nan(2,n_array); % store the # of indicies above and below the center value
    for jj = 1:n_array % for each column
        [cIdx(1),~] = find( data(:,jj)==cVal(jj) ); % find location of center value
        R(1,jj) = length(1:cIdx(1)-1); % length of data below ceneter
        R(2,jj) = length((cIdx(1)+1):n_center); % length of data above center
    end

    dR = diff(R); % change in length needed to center data for each column

    n_new = n_center + max(abs(dR)); % new length
    centData = nan(n_new,n_array); % centered data
    for jj = 1:n_array % for each column
        if dR(jj)>=1
            centData(:,jj) = [nan(dR(jj),1) ; data(:,jj)];
        elseif dR(jj)<=-1
            centData(:,jj) = [data(:,jj) ; nan(abs(dR(jj)),1)];
        elseif dR(jj)==0
            centData(:,jj) = data(:,jj);
        end
    end

    if dim==2 % transform back to row vectors
        R = R';
        dR = dR';
        centData = centData';
    end
else
    Error('Error: input must be matrix or cell array') 
end
end