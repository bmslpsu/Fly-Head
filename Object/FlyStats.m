classdef FlyStats
    %% Stats:
    % 
    properties (GetAccess=private) % properties only avaiable to class
    end
    
    properties (Constant) % constant properties
    end
   
    properties (Dependent) % dependent properties
    end
    
  	properties % all other properties
        All         = [];       % all objects
        Mean        = [];       % mean of objects
        Median   	= [];       % median of objects
        STD         = [];       % std of objects
        Var         = [];       % variance of objects
        Max         = [];       % max of objects
        Min         = [];       % min of objects
        Mode        = [];       % mode of objects
        Range       = [];       % range of objects
        CircMean  	= [];       % circular mean of objects
        CircSTD   	= [];       % circular std of objects
    end
    
    methods
        function obj = FlyStats(obj_cell)
            [nn] = size(obj_cell); % size of cell array

            OBJ = cell(nn(1),nn(2)); % stores object data
            % Transform objects into cell arrays
            for kk = 1:nn(1) % cycle through all input cells in the dimesnion to average
                for jj = 1:nn(2) % cycle through all input cells in other dimension
                prop = properties(obj_cell{kk,jj}); % object property names
                    for ii = 1:length(prop) % cycle through properties
                      name = prop{ii}; % get property name
                      OBJ{kk,jj}{ii,1} = obj_cell{kk,jj}.(name); % get value associated with property & store in cell array
                    end
                end
            end

            obj_all = cell(1,nn(2)); % stores concatenateed cell data
            obj.All = cell(1,nn(2)); % stores concatenateed array data
            for jj = 1:nn(2) % cycle through all input cells in the dimension average
                obj_all{jj} = cat(3,OBJ{:,jj}); % concatenates object data cells in 3rd dimension
                for ii = 1:length(obj_all{jj}) % cycle through properties
                    for kk = 1:nn(1) % cycle through the dimesnion to average
                        obj.All{ii,jj}(:,:,kk) = obj_all{jj}{ii,1,kk}; % transform cell to array
                    end
                end
            end
            
            % Statistics
            obj.Mean        = cellfun(@(x) mean(x,3),           obj.All,'UniformOutput',false);
            obj.Median      = cellfun(@(x) median(x,3),         obj.All,'UniformOutput',false);
            obj.STD         = cellfun(@(x) std(x,0,3),          obj.All,'UniformOutput',false);
            obj.Var         = cellfun(@(x) var(x,0,3),          obj.All,'UniformOutput',false);
            obj.Max         = cellfun(@(x) max(x,[],3),         obj.All,'UniformOutput',false);
            obj.Min         = cellfun(@(x) min(x,[],3),         obj.All,'UniformOutput',false);
            obj.Mode        = cellfun(@(x) mode(x,3),           obj.All,'UniformOutput',false);
            obj.Range       = cellfun(@(x) range(x,3),          obj.All,'UniformOutput',false);
        	obj.CircMean    = cellfun(@(x) circ_mean(x,[],3),   obj.All,'UniformOutput',false);
            obj.CircSTD     = cellfun(@(x) circ_std(x,[],[],3), obj.All,'UniformOutput',false);

        end
    end
end