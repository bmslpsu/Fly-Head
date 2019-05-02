classdef GrandStats
    %% GrandStats: statistics for multiple instances of "Fly" or "IO" classes
    %   INPUTS:
    %       FlyStat    : cell array containing FlyStat classes (works along columns)
    properties (GetAccess=private) % properties only avaiable to class
        nFly   	= [];           % # of flies
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
            obj.nFly = nn(1);
            
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
        
        function [] = PlotTime(obj,n)
            % PlotTime: time domain plot
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            if nargin==1
                n = 1:size(obj.All{1},2);
            end
            
            figure('Name',['Grand Time Domain (# Fly = ' num2str(obj.nFly) ')'])
            pp = 1;
            time = obj.Median{2}{5};
            for kk = n
                ax = subplot(length(n),1,pp) ; hold on ; grid on
                ylabel(['X_' num2str(kk)])
                
                if kk==n(end)
                    xlabel('Time')
                end
                xlim([0 time(end)])
                
                for jj = 1:size(obj.All{2}{5},3)
                    h = plot(obj.All{2}{5}(:,1,jj),obj.All{2}{6}(:,kk,jj));
                    h.Color(4) = 0.5;
                end
                
                PlotPatch(obj.Median{2}{6}(:,kk),obj.STD{2}{6}(:,kk),time,2,obj.nFly,'k',[0.4 0.4 0.6],0.5,2);
                
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
        
            freq = obj.Median{2}{7}(:,1);
            if nargin==1
                n = 1:size(obj.All{6},2);
                lim = freq(end);
            elseif nargin<=2
                lim = freq(end);
            end
            
            figure('Name',['Grand Frequency Domain (# Fly = ' num2str(obj.nFly) ')'])
            pp = 1;
            for kk = n
                ax = subplot(2,length(n),pp) ; hold on ; grid on
                    title(['X_' num2str(kk)])               
                    for jj = 1:size(obj.All{1}{1},3)
                        h = plot(obj.All{2}{7}(:,1,jj),obj.All{2}{8}(:,kk,jj));
                        h.Color(4) = 0.5;
                       	h = plot(obj.All{2}{10}(:,1,jj),obj.All{2}{11}(:,kk,jj));
                        h.Color(4) = 0.5;
                    end

                    PlotPatch(obj.Median{2}{8}(:,kk),obj.STD{2}{8}(:,kk),freq,2,obj.nFly,'k',[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.STD{2}{11}(:,kk);
                    errorbar(obj.Median{2}{10},obj.Median{2}{11}(:,kk),err,'-or','LineWidth',2)

                    xlim([0 lim])
                    ylim(max(ax.YLim)*[0 1])
                    ylabel('Magnitude')
                
                ax = subplot(2,length(n),pp+length(n)) ; hold on ; grid on                               
                    for jj = 1:size(obj.All{5},3)
                        h = plot(obj.All{2}{7}(:,1,jj),obj.All{9}{9}(:,kk,jj));
                        h.Color(4) = 0.5;
                     	h = plot(obj.All{2}{10}(:,1,jj),obj.All{9}{12}(:,kk,jj));
                        h.Color(4) = 0.5;
                    end

                    PlotPatch(obj.CircMean{9}{9}(:,kk),obj.CircSTD{9}{9}(:,kk),freq,2,obj.nFly,'k',[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.CircSTD{9}{12}(:,kk);
                    errorbar(obj.Median{2}{10},obj.CircMean{9}{12}(:,kk),err,'-or','LineWidth',2)
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
            % PlotIOBode: BODE plot
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
        
            freq = obj.Median{2}{7}(:,1);
            if nargin==1
                n = 1:size(obj.All{6},2);
                lim = freq(end);
            elseif nargin<=2
                lim = freq(end);
            end
            
            figure('Name',['Grand BODE (# Fly = ' num2str(obj.nFly) ')'])
            pp = 1;
            for kk = n
                ax = subplot(2,length(n),pp) ; hold on ; grid on
                    title(['X_' num2str(kk)])               
                    for jj = 1:size(obj.All{1}{1},3)
                       	h = plot(obj.All{2}{4}(:,1,jj),obj.All{2}{5}(:,kk,jj),'-o');
                        h.Color(4) = 0.5;
                    end

                    err = 2*obj.STD{2}{5}(:,kk);
                    errorbar(obj.Median{2}{4},obj.Median{2}{5}(:,kk),err,'-ok','LineWidth',2)

                    xlim([0 lim])
                    ylim(max(ax.YLim)*[0 1])
                    ylabel('Gain')
                
                ax = subplot(2,length(n),pp+length(n)) ; hold on ; grid on                               
                    for jj = 1:size(obj.All{1}{1},3)
                     	h = plot(obj.All{2}{4}(:,1,jj),obj.All{9}{6}(:,kk,jj),'-o');
                        h.Color(4) = 0.5;
                    end

                    err = 2*obj.CircSTD{9}{6}(:,kk);
                    errorbar(obj.Median{2}{4},obj.CircMean{9}{6}(:,kk),err,'-ok','LineWidth',2)
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
        
            freq = obj.Median{2}{8}(:,1);
            if nargin==1
                n = 1:size(obj.All{2}{6},2);
                lim = freq(end);
            elseif nargin<=2
                lim = freq(end);
            end
            
            figure('Name',['Grand Coherence (# trial = ' num2str(obj.nFly) ')'])
            pp = 1;
            for kk = n
                subplot(1,length(n),pp) ; hold on ; grid on
                    title(['X_' num2str(kk)])               
                    for jj = 1:size(obj.All{2}{1},3)
                        h = plot(obj.All{2}{8}(:,1,jj),obj.All{2}{7}(:,kk,jj),'-');
                        h.Color(4) = 0.5;
                        cc = h.Color;
                        h = plot(obj.All{2}{4}(:,1,jj),obj.All{2}{9}(:,kk,jj),'-o');
                        h.Color(:) = cc;
                    end
                    
                    PlotPatch(obj.Median{2}{7}(:,kk),obj.STD{2}{7}(:,kk),freq,2,obj.nFly,[0.5 0.5 0.5],[0.4 0.4 0.6],0.5,2);

                    err = 2*obj.STD{2}{9}(:,kk);
                    errorbar(obj.Median{2}{4},obj.Median{2}{9}(:,kk),err,'-ok','LineWidth',2)

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