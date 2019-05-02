classdef FlyStats
    %% FlyStats: statistics for multiple instances of "Fly" or "IO" classes
    %   INPUTS:
    %       obj_cell    : cell array containing classes (works along columns so columns must contain the same class)
    properties (GetAccess=private) % properties only avaiable to class
        nTrial    	= [];       % # of trials
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
            obj.nTrial = nn(1);

            OBJ = cell(nn(1),nn(2)); % stores object data
            % Transform objects into cell arrays
            for kk = 1:nn(1) % cycle through all input cells in the dimension to average
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
        
        function nTrial = get.nTrial(obj)
            nTrial = obj.nTrial;
         end
        
        
        function [] = PlotTime(obj,n)
            % PlotTime: time domain plot
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            if nargin==1
                n = 1:size(obj.All{6},2);
            end
            
            figure('Name',['Fly Time Domain (# trial = ' num2str(obj.nTrial) ')'])
            pp = 1;
            time = obj.Median{5};
            for kk = n
                ax = subplot(length(n),1,pp) ; hold on ; grid on
                ylabel(['X_' num2str(kk)])
                
                if kk==n(end)
                    xlabel('Time')
                end
                xlim([0 time(end)])
                
                for jj = 1:size(obj.All{5},3)
                    h = plot(obj.All{5}(:,1,jj),obj.All{6}(:,kk,jj));
                    h.Color(4) = 0.5;
                end
                
                PlotPatch(obj.Median{6}(:,kk),obj.STD{6}(:,kk),time,2,obj.nTrial,'k',[0.4 0.4 0.6],0.5,2);
                
                ylim(max(ax.YLim)*[-1 1])
                
                pp = pp + 1;
            end            
            hold off
        end
        
    	function [] = PlotFreq(obj,n,lim)
            % PlotFreq: frequency domain plot
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
        
            freq = obj.Median{7}(:,1);
            if nargin==1
                n = 1:size(obj.All{6},2);
                lim = freq(end);
            elseif nargin<=2
                lim = freq(end);
            end
            
            figure('Name',['Fly Frequency Domain (# trial = ' num2str(obj.nTrial) ')'])
            pp = 1;
            for kk = n
                ax = subplot(2,length(n),pp) ; hold on ; grid on
                    title(['X_' num2str(kk)])               
                    for jj = 1:size(obj.All{5},3)
                        h = plot(obj.All{7}(:,1,jj),obj.All{8}(:,kk,jj));
                        h.Color(4) = 0.5;
                        h = plot(obj.All{10}(:,1),obj.All{11}(:,kk,jj),'-o');
                        h.Color(4) = 0.5;
                    end

                    PlotPatch(obj.Median{8}(:,kk),obj.STD{8}(:,kk),freq,2,obj.nTrial,'k',[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.STD{11}(:,kk);
                    errorbar(obj.Median{10},obj.Median{11}(:,kk),err,'-or','LineWidth',2)

                    xlim([0 lim])
                    ylim(max(ax.YLim)*[0 1])
                    ylabel('Magnitude')
                
                ax = subplot(2,length(n),pp+length(n)) ; hold on ; grid on                               
                    for jj = 1:size(obj.All{5},3)
                        h = plot(obj.All{7}(:,1,jj),obj.All{9}(:,kk,jj));
                        h.Color(4) = 0.5;
                     	h = plot(obj.All{10}(:,1),obj.All{12}(:,kk,jj),'-o');
                        h.Color(4) = 0.5;
                    end

                    PlotPatch(obj.CircMean{9}(:,kk),obj.CircSTD{9}(:,kk),freq,2,obj.nTrial,'k',[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.CircSTD{12}(:,kk);
                    errorbar(obj.Median{10},obj.CircMean{12}(:,kk),err,'-or','LineWidth',2)

                    xlim([0 lim])
                    ylim(max(ax.YLim)*[-1 1])
                    ylabel('Phase')

                    if (pp+length(n))>n(end)
                        xlabel('Frequency')
                    end
                
                pp = pp + 1;
            end
            hold off
        end
        
        
        function [] = PlotIOBode(obj,n,lim)
            % PlotIOBode: BODE plot for IO_Class
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
        
            freq = obj.Median{1}(:,1);
            if nargin==1
                n = 1:size(obj.All{6},2);
                lim = freq(end);
            elseif nargin<=2
                lim = freq(end);
            end
            
            figure('Name',['Fly Frequency Domain (# trial = ' num2str(obj.nTrial) ')'])
            pp = 1;
            for kk = n
                ax = subplot(2,length(n),pp) ; hold on ; grid on
                    title(['X_' num2str(kk)])               
                    for jj = 1:size(obj.All{1},3)
                        h = plot(obj.All{4}(:,1,jj),obj.All{5}(:,kk,jj),'-o');
                        h.Color(4) = 0.5;
                    end

%                     PlotPatch(obj.Median{2}(:,kk),obj.STD{2}(:,kk),freq,2,obj.nTrial,'k',[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.STD{5}(:,kk);
                    errorbar(obj.Median{4},obj.Median{5}(:,kk),err,'-ok','LineWidth',2)

                    xlim([0 lim])
                    ylim([0 1.2])
                    ylabel('Gain')
                
                ax = subplot(2,length(n),pp+length(n)) ; hold on ; grid on                               
                    for jj = 1:size(obj.All{5},3)
                        h = plot(obj.All{4}(:,1,jj),obj.All{6}(:,kk,jj),'-o');
                        h.Color(4) = 0.5;
                    end

                    err = 2*obj.CircSTD{6}(:,kk);
                    errorbar(obj.Median{4},obj.CircMean{6}(:,kk),err,'-ok','LineWidth',2)

                    xlim([0 lim])
                    ylim(max(ax.YLim)*[-1 1])
                    ylabel('Phase')

                    if (pp+length(n))>n(end)
                        xlabel('Frequency')
                    end
                
                pp = pp + 1;
            end
            hold off
        end
        
    function [] = PlotCohr(obj,n,lim)
            % PlotCohr: Coherence plot for IO_Class
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
        
            freq = obj.Median{8}(:,1);
            if nargin==1
                n = 1:size(obj.All{6},2);
                lim = freq(end);
            elseif nargin<=2
                lim = freq(end);
            end
            
            figure('Name',['Fly Coherence (# trial = ' num2str(obj.nTrial) ')'])
            pp = 1;
            for kk = n
                subplot(1,length(n),pp) ; hold on ; grid on
                    title(['X_' num2str(kk)])               
                    for jj = 1:size(obj.All{1},3)
                        h = plot(obj.All{8}(:,1,jj),obj.All{7}(:,kk,jj),'-');
                        h.Color(4) = 0.5;
                        cc = h.Color;
                        h = plot(obj.All{4}(:,1,jj),obj.All{9}(:,kk,jj),'-o');
                        h.Color(:) = cc;
                    end
                    
                    PlotPatch(obj.Median{7}(:,kk),obj.STD{7}(:,kk),freq,2,obj.nTrial,[0.5 0.5 0.5],[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.STD{9}(:,kk);
                    errorbar(obj.Median{4},obj.Median{9}(:,kk),err,'-ok','LineWidth',2)

                    xlim([0 lim])
                    ylim([0 1])
                    ylabel('Coherence')
                    
                    if pp==n(end)
                        xlabel('Frequency')
                    end
                
                pp = pp + 1;
            end
            hold off
        end
        
    end
end