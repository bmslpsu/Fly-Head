function [D,I,N,U] = GetFileData(FILES,varargin)
%% GetFileData: Parse file name data and returns tables with relevant information, lets user load files if no input is specified
%   INPUTS:
%       FILES       :   file cells in the form "var1_val1_var2_val2_..._varn_valn". The first variable is
%                       the control group & the second is the repetions, the rest are categories
%       varargin    :   user can rename the variables by inputing strings equal to the # of variables in
%                       the file name
%   OUTPUTS:
%       D           :   raw file data table
%       I           :   index table
%       N           :   # table
%       U           :   unique variable #'s table
%---------------------------------------------------------------------------------------------------------------------------------
% Let user load files if no input is specified
if ~nargin
    [files, ~] = uigetfile({'*', 'files'}, 'Select files', 'MultiSelect','on');
    FILES = cellstr(files)';
end
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
% Compare whther data is string or number >> determine if category or value
log = nan(1,n.vars);
for ii = 1:n.vars
    log(ii) = isnan(str2double(vardata{1,ii}));
end
if ~log(1) % make sure naming convetion is consistent and file name starts with a category
    warning('file name should not start with a number')
end
for kk = 1:length(log)-1
   if (log(kk)+log(kk+1))==0 % make sure there are not values without categories
      warning('file naming-convention is not correct')
   end
end
[~,loc.catg] = find(log==true); % find locations of categorical variable 
[~,loc.val] = find(log==false); % find locations of value variables
loc.catg = loc.catg(1:length(loc.val)); % if there are more categories than values >>> get rid of last category
n.catg = length(loc.val); % # of categories
n.val = length(loc.catg); % # of values
if length(varargin)>n.catg
   error('More variable names than file variables') 
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
    vardata(:,loc.val(kk)) = cellfun(@(x) str2double(x),vardata(:,loc.val(kk)),...
        'UniformOutput',false); % get values from cell array
    numdata(:,kk) = cell2mat(vardata(:,loc.val(kk))); % numeric array of cvategory values
    unq{kk} = sort(unique(numdata(:,kk)),'ascend'); % unique values for category
    nn(kk) = length(unq{kk}); % # of values for category
    idx{kk} = (1:nn(kk))'; % indicies values needed for each category
    % Set up indexing convention
    reps{kk} = nan(nn(kk),1);
    for jj = 1:nn(kk)
        reps{kk}(jj) = sum(numdata(:,kk)==unq{kk}(jj)); % # of repetitions of each unique value
        Ind(numdata(:,kk)==unq{kk}(jj),kk) = jj; % set idicies to start at 1 & increment by 1
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

% Let user set variable names if needed
varnames = catg;
if nargin>1
    for kk = 1:nargin-1
        varnames{kk} = varargin{kk};
    end
end

% Find map for conditions
cond = nn(3:end); % conditions per category
ncatg = length(cond); % # of categories 
ridx = cell(ncatg,1);
for kk = 1:ncatg
    ridx{kk} = idx{2+kk}';
end
allcomb = combvec(ridx{:}); % all unquie condition combinations
map = nan(nn(1),size(allcomb,2)); % map to store # reps for each combination
for kk = 1:nn(1) % per fly
    % Check all conditions
    for ii = 1:size(allcomb,2) % each condition
        rr = ones(size(Ind,1),1);
        for jj = 1:size(allcomb,1) % each condition
            rr = rr & (Ind(:,2+jj)==allcomb(jj,ii));
        end
        map(kk,ii) = sum( rr & Ind(:,1)==kk ); % # of reps
    end
end

% Check if there are at least 3 reps per condition
minTrial = 3;
mpass = nan(nn(1),1);
for kk = 1:nn(1)
    mcheck = min(map(kk,:));
    if mcheck>=minTrial
        mpass(kk) = true;
    else
        mpass(kk) = false;
    end
end
npass = sum(mpass); % # of successful flies
mpass = logical(mpass);

% Make variable names for map
mapname = cell(1,size(allcomb,2));
for kk = 1:size(allcomb,2)
    valstr = [];
    for jj = 1:size(allcomb,1)
       valstr = [valstr  '_' num2str(allcomb(jj,kk))];
    end
    mapname{kk} = [strcat(varnames{3:end})  valstr];
end

% Make table from raw file data
D = splitvars(table(numdata));
D.Properties.VariableNames = varnames;
% Make table from index data
I = splitvars(table(Ind));
I.Properties.VariableNames = varnames;
% Make table from # data
N = splitvars(table([nn,n.file]));
N.Properties.VariableNames = [varnames,'file'];
% Make table from unqiue data
U = splitvars(table(unq));
U.Properties.VariableNames = varnames;

fprintf('Files: %i \n',n.file)
fprintf('%s: %i \n',varnames{1},nn(1))
for kk = 3:n.catg
    fprintf('%s: %i \n',varnames{kk},nn(kk))
end
fprintf('Pass: %i \n',npass)

if isempty(map)
    mpass = reps{:,1} >= minTrial;
    T = splitvars(table(idx{1} , reps{:,1} , mpass));
    T.Properties.VariableNames = [varnames(1:2) , ['CHECK_' num2str(minTrial)]];
else
    T = splitvars(table(idx{1} , reps{:,1} , map , mpass));
    T.Properties.VariableNames = [varnames(1:2) , mapname ,['CHECK_' num2str(minTrial)]];
end
disp(T)
end