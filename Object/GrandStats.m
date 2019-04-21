classdef GrandStats
    %% GrandStats:
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
        function obj = GrandStats(FlyStat)
            [nn] = size(FlyStat); % size of cell array
            
%             obj.All = cell(nn(1),nn(2));
            for ii = 1:length(FlyStat{1}.Mean)
                for kk = 1:nn(1)
                    prop = properties(FlyStat{kk});
                    for qq = 2:length(prop)
                        name = prop{qq}; % get property name
                        for ww = 1:size(FlyStat{kk}.(name),2)
                            obj.All{qq-1,1}{ii,ww}(:,:,kk) = FlyStat{kk}.(name){ii,ww};
                        end
                    end
                end
            end
            
            % Statistics
            for ii = 1:length(obj.All)
                obj.Mean        {ii,1}   	= cellfun(@(x) mean(x,3),           obj.All{ii},'UniformOutput',false);
                obj.Median      {ii,1}      = cellfun(@(x) median(x,3),         obj.All{ii},'UniformOutput',false);
                obj.STD         {ii,1}   	= cellfun(@(x) std(x,0,3),          obj.All{ii},'UniformOutput',false);
                obj.Var         {ii,1}   	= cellfun(@(x) var(x,0,3),          obj.All{ii},'UniformOutput',false);
                obj.Max         {ii,1}      = cellfun(@(x) max(x,[],3),         obj.All{ii},'UniformOutput',false);
                obj.Min         {ii,1}  	= cellfun(@(x) min(x,[],3),         obj.All{ii},'UniformOutput',false);
                obj.Mode        {ii,1}  	= cellfun(@(x) mode(x,3),           obj.All{ii},'UniformOutput',false);
                obj.Range       {ii,1}      = cellfun(@(x) range(x,3),          obj.All{ii},'UniformOutput',false);
                obj.CircMean 	{ii,1}   	= cellfun(@(x) circ_mean(x,[],3),  	obj.All{ii},'UniformOutput',false);
                obj.CircSTD 	{ii,1}   	= cellfun(@(x) circ_std(x,[],[],3),	obj.All{ii},'UniformOutput',false);
            end      


        end
    end
end