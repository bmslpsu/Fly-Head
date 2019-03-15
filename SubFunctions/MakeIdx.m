function [FD] = MakeIdx(FILES,varargin)
%% MakeIdx: 
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
%---------------------------------------------------------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------------------------------------------------
n.file = length(FILES);
for kk = 1:n.file
    [~,filename,~] = fileparts(FILES{kk}); % remove file extension from file names
    filedata = textscan(char(filename), '%s', 'delimiter', '_');
    filedata = filedata{1}';
    n.vars = length(filedata);
%     vardata = cell(n.file,n.vars);
    for jj = 1:length(filedata)
        vardata{kk,jj} = filedata{jj};
    end
end
log = nan(1,n.vars);
for ii = 1:n.vars
    log(ii) = isnan(str2double(vardata{1,ii}));
end
[~,var] = find(log==false);
[~,str] = find(log==true);
str = str(1:length(var));
name = cell(1,length(str));
% num  = cell(1,length(var));
for ii = 1:length(str)
    name{ii} = vardata{1,str(ii)};
% 	num {ii} = vardata{1,var(ii)};
end

for kk = 1:length(var)
    unq(kk) = unique(cell2mat(vardata{:,var(kk)}));
end




end