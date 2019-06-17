function [D,I,N,U,T,FILES,PATH] = GetFileData(FILES,abscat,varargin)
%% GetFileData: Parse file name data and returns tables with relevant information
%   INPUTS:
%       FILES       :   file cells in the form "var1_val1_var2_val2_..._varn_valn". The first variable is
%                       the control group & the second is the repetions,
%                       the rest are categories. If "FILES" is a character
%                       vector, then it is the root directory & user selects files.
%       abscat      :   if set to true, indexing is based on the absolute
%                       value of conditions.
%       varargin    :   user can rename the variables by inputing strings equal to the # of variables in
%                       the file name.
%   OUTPUTS:
%       D           :   raw file data table
%       I           :   index table
%       N           :   # table
%       U           :   unique variable #'s table
%       T           :   map of all data + trials per condition
%       FILES      	:   files used
%       PATH      	:   file location
%---------------------------------------------------------------------------------------------------------------------------------
% Let user load files if no input is specified
if ~nargin % open file selection GUi in current folder
    [files, PATH] = uigetfile({'*', 'files'}, 'Select files', 'MultiSelect','on');
    FILES = cellstr(files)';
    abscat = false;
elseif nargin>=1
    if ischar(FILES) % if root directory given, open file selection GUi in root
        [files, PATH] = uigetfile({'*', 'files'}, 'Select files',FILES, 'MultiSelect','on');
        FILES = cellstr(files)';
    else
        PATH = []; % if only files are given
    end
    
    if nargin==1
        abscat = false; % default is off
    end
end
clear files

% Get file data
n.file = length(FILES);
[~,filename,~] = fileparts(FILES{1}); % remove file extension from file name
filedata = textscan(char(filename), '%s', 'delimiter', '_'); % deliminate file name
n.vars = length(filedata{1}); % # of variables in file name
vardata = cell(n.file,n.vars);
for kk = 1:n.file
    [~,filename,~] = fileparts(FILES{kk}); % remove file extension from file name
    filedata = textscan(char(filename), '%s', 'delimiter', '_'); % deliminate file name
    filedata = filedata{1}'; % get data from cell
    n.vars = length(filedata); % # of variables in file name
    % Store file data in cell array
    for jj = 1:length(filedata)
        vardata{kk,jj} = filedata{jj};
    end
end

% Compare whether data is string or number >> determine if category or value
log = false(1,n.vars);
for ii = 1:n.vars
    log(ii) = isnan(str2double(vardata{1,ii}));
end

if ~log(1) % make sure naming convetion is consistent and file name starts with a category
    warning('file name should not start with a number')
end

for kk = 1:length(log)-1
   if (log(kk)+log(kk+1))==0 % make sure there are not values without categories
      warning('file naming convention is not correct')
   end
   
   if log(kk+1)==1 % if we end on a category, use it as an index
        catFlag = true;
   else
        catFlag = false;
   end
end

[~,loc.catg] = find(log==true); % find locations of categorical variable 
[~,loc.val] = find(log==false); % find locations of value variables
loc.catg = loc.catg(1:length(loc.val)); % if there are more categories than values >>> get rid of trailing categories
n.catg = length(loc.val); % # of categories
n.val = length(loc.catg); % # of values
if length(varargin)>n.catg
   error('More variable names than variables')
end

catg = cell(1,length(loc.catg)); % cell to store category names
for ii = 1:n.val
    catg{ii} = vardata{1,loc.catg(ii)}; % get names
end

% Store values, get unique values & # of unique values
numdata = nan(n.file,n.val);
unq = cell(1,n.val);
nn = nan(1,n.val);
idx = cell(1,n.val);
reps = cell(1,n.val);
Ind = nan(n.file,n.val);
for kk = 1:n.val
    if abscat
        vardata(:,loc.val(kk)) = cellfun(@(x) abs(str2double(x)),vardata(:,loc.val(kk)),...
            'UniformOutput',false); % get absolute values from cell array
    else
        vardata(:,loc.val(kk)) = cellfun(@(x) str2double(x),vardata(:,loc.val(kk)),...
            'UniformOutput',false); % get values from cell array
    end
    
    numdata(:,kk) = cell2mat(vardata(:,loc.val(kk))); % numeric array of category values
    unq{kk} = abs_sort(unique(numdata(:,kk))); % unique values for category
    
    try
        nanIdx = find(isnan(unq{3}));
    catch
        nanIdx = [];
    end
    
    if ~isempty(nanIdx) % make sure NaN only registers once
        unq{kk}(nanIdx(2:end)) = [];
    end
    
    nn(kk) = length(unq{kk}); % # of values for category
    idx{kk} = (1:nn(kk))'; % indicie values needed for each category
    
    % Set up indexing convention
    reps{kk} = nan(nn(kk),1);
    for jj = 1:nn(kk)
        reps{kk}(jj) = sum(numdata(:,kk)==unq{kk}(jj)); % # of repetitions of each unique value
        Ind(numdata(:,kk)==unq{kk}(jj),kk) = jj; % set idicies to start at 1 & increment by 1
        if isnan(unq{kk}(jj)) % replace Nan with index
            [rr,~] = find(isnan(numdata(:,kk)));
            Ind(rr,kk) = jj;
        end
    end
end

% Associate trials with main index (first category)
pp = 1;
tt = 2; % trial index
for jj = 1:nn(1)
    rng = pp:pp+reps{1}(jj)-1; % range of trials that main index spans
    Ind(rng,tt) = 1:reps{1}(jj); % set idicies to start at 1 & increment by 1
    pp = pp + reps{1}(jj); % new start of range
end

% If filenames ends on a category, use it as an index
if catFlag
  	catIdx  = nan(n.file,1);
   	catData = vardata(:,end);
    catUnq  = sort(unique(catData));
    catN    = length(catUnq);
    for kk = 1:catN
        catLoc = strcmp(catData,catUnq{kk});
      	catIdx(catLoc) = kk;
    end
	numdata(:,end+1)	= catIdx;
    Ind(:,end+1)        = catIdx;
    nn(end+1)           = catN;
	idx{end+1}          = (1:catN)';
    unq{end+1}          = catUnq;
    catg{end+1}         = 'class';
end

% Let user set variable names if specififed
varnames = catg;
if nargin>2
    for kk = 1:nargin-2
        varnames{kk} = varargin{kk};
    end
end

% Compute trial map for conditions
cond  = nn(3:end); % conditions per category
ncatg = length(cond); % # of categories 
ridx  = cell(ncatg,1);
for kk = 1:ncatg
    ridx{kk} = idx{2+kk}';
end
allcomb = combvec(ridx{:}); % all unique condition combinations
ncomb = size(allcomb,2);
map = nan(nn(1),ncomb); % map to store # reps for each combination
for kk = 1:nn(1) % per fly
    % Check all combinations & conditions
    for ii = 1:ncomb % each combination
        rr = ones(size(Ind,1),1);
        for jj = 1:ncatg % each condition
            rr = rr & (Ind(:,2+jj)==allcomb(jj,ii));
        end
        map(kk,ii) = sum( rr & Ind(:,1)==kk ); % # of reps
    end
end

% If there are no condtions, use trials instead >>>  make variable names for map
if isempty(map) % no conditions
   map = reps{:,1};
   mapname = {'reps'};
else
    % Make map variable names
    mapname = cell(1,ncomb);
    for kk = 1:ncomb
        valstr = cell(1,ncomb);
        for jj = 1:ncatg
            if jj==ncatg
                valstr{jj} = strcat(varnames{2+jj},['_' num2str(allcomb(jj,kk))]);
            else
                valstr{jj} = strcat(varnames{2+jj},['_' num2str(allcomb(jj,kk)) '_']);
            end
        end
        mapname{kk} = cat(2,valstr{:});
    end
end

% Check if there are at least 3 reps per condition & if each fly has at least one rep
minTrial = 3;
mpass = nan(nn(1),1);
for kk = 1:nn(1)
    mcheck = min(map(kk,:));
    if mcheck>=minTrial
        mpass(kk) = true;
    else
        mpass(kk) = false;
    end
    
    if ~mcheck % makes sure at least one trial per condition
        warning('Fly %i does not have data for all conditions',unq{1}(kk))
    end
end
npass = sum(mpass); % # of successful flies
mpass = logical(mpass);
mapSum = sum(map,1); % # of trials per condition

% Table from raw file data
D = splitvars(table(numdata));
D.Properties.VariableNames = varnames;
% Table from index data
I = splitvars(table(Ind));
I.Properties.VariableNames = varnames;
% Table from # data
N = splitvars(table([nn,n.file]));
N.Properties.VariableNames = [varnames,'file'];
% Ttable from unique data
U = splitvars(table(unq));
U.Properties.VariableNames = varnames;
% Table from map conditions
M = splitvars(table(mapSum));
M.Properties.VariableNames = mapname;
% Table from all trials ,flies, conditions
T = splitvars(table(unq{1} , reps{:,1} , map , mpass));
T.Properties.VariableNames = [varnames(1:2) , mapname ,['CHECK_' num2str(minTrial)]];

% Display data
fprintf('Files: %i \n',n.file)
fprintf('%s: %i \n',varnames{1},nn(1))
for kk = 3:n.catg
    fprintf('%s: %i \n',varnames{kk},nn(kk))
end
fprintf('Pass: %i \n',npass)
disp(T)
disp(M)
end