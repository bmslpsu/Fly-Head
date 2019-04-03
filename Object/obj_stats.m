function [data] = obj_stats(obj_cell)
%% obj_stats: computes statistics of objects of the same class stored in a cell array
% Only works along 1st dimension for now
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------
[nn] = size(obj_cell); % size of cell array

obj = cell(nn(1),nn(2)); % stores object data
% Transform objects into cell arrays
for kk = 1:nn(1) % cycle through all input cells in the dimesnion to average
    for jj = 1:nn(2) % cycle through all input cells in other dimension
    prop = properties(obj_cell{kk,jj}); % object property names
        for ii = 1:length(prop) % cycle through properties
          name = prop{ii}; % get property name
          obj{kk,jj}{ii,1} = obj_cell{kk,jj}.(name); % get value associated with property & store in cell array
        end
    end
end

obj_all = cell(1,nn(2)); % stores concatenateed cell data
data.all = cell(1,nn(2)); % stores concatenateed array data
for jj = 1:nn(2) % cycle through all input cells in the dimesnion average
	obj_all{jj} = cat(3,obj{:,jj}); % concatenates object data cells in 3rd dimension
    for ii = 1:length(obj_all{jj}) % cycle through properties
        for kk = 1:nn(1) % cycle through the dimesnion to average
            data.all{jj}{ii}(:,:,kk) = obj_all{jj}{ii}; % transform cell to array
        end
    end
end

% Statistics calculations
data.mean   = cell(1,nn(2));
data.median = cell(1,nn(2));
data.std    = cell(1,nn(2));
data.var    = cell(1,nn(2));
data.max    = cell(1,nn(2));
data.min    = cell(1,nn(2));
data.mode 	= cell(1,nn(2));
data.range  = cell(1,nn(2));
for jj = 1:nn(2)
    for ii = 1:length(data.all{jj})
        data.mean       {jj}{ii} = mean     (data.all{jj}{ii},3);
        data.median     {jj}{ii} = median   (data.all{jj}{ii},3);
        data.std        {jj}{ii} = std      (data.all{jj}{ii},0,3);
        data.var        {jj}{ii} = var      (data.all{jj}{ii},0,3);
        data.max        {jj}{ii} = max      (data.all{jj}{ii},[],3);
        data.min        {jj}{ii} = min      (data.all{jj}{ii},[],3);
        data.mode     	{jj}{ii} = mode     (data.all{jj}{ii},3);
        data.range     	{jj}{ii} = range  	(data.all{jj}{ii},3);
    end
end

end