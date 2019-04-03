classdef Stats
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
    end
    
    methods
        function obj = Stats(obj_cell)
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
            for jj = 1:nn(2) % cycle through all input cells in the dimesnion average
                obj_all{jj} = cat(3,OBJ{:,jj}); % concatenates object data cells in 3rd dimension
                for ii = 1:length(obj_all{jj}) % cycle through properties
                    for kk = 1:nn(1) % cycle through the dimesnion to average
                        obj.All{jj}{ii}(:,:,kk) = obj_all{jj}{ii}; % transform cell to array
                    end
                end
            end            

            % Statistics calculations
            obj.Mean   = cell(1,nn(2));
            obj.Median = cell(1,nn(2));
            obj.STD    = cell(1,nn(2));
            obj.Var    = cell(1,nn(2));
            obj.Max    = cell(1,nn(2));
            obj.Min    = cell(1,nn(2));
            obj.Mode 	= cell(1,nn(2));
            obj.Range  = cell(1,nn(2));
            for jj = 1:nn(2)
                for ii = 1:length(obj.All{jj})
                    obj.Mean        {jj}{ii} = mean     (obj.All{jj}{ii},3);
                    obj.Median      {jj}{ii} = median   (obj.All{jj}{ii},3);
                    obj.STD         {jj}{ii} = std      (obj.All{jj}{ii},0,3);
                    obj.Var         {jj}{ii} = var      (obj.All{jj}{ii},0,3);
                    obj.Max         {jj}{ii} = max      (obj.All{jj}{ii},[],3);
                    obj.Min         {jj}{ii} = min      (obj.All{jj}{ii},[],3);
                    obj.Mode     	{jj}{ii} = mode     (obj.All{jj}{ii},3);
                    obj.Range   	{jj}{ii} = range  	(obj.All{jj}{ii},3);
                end
            end
        end
    end
end