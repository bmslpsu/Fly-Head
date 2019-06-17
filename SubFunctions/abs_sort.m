function [sortData] = abs_sort(data)
%% abs_sort: sort a vector by seperating the sorted positive & negative components and concatenating the results
%   INPUTS:
%       data        :   data to sort
%   OUTPUTS:
%       sortData	:   sorted data
%---------------------------------------------------------------------------------------------------------------------------------
dim = size(data);
dIdx = dim==1;
if sum(dIdx)>1 ||sum(dIdx)==0
    error('Error: input must be vector')
end
dim = find(dIdx==0);

pIdx = data>=0;
nIdx = data<0;

pSort = sort(data(pIdx),'ascend');
nSort = sort(data(nIdx),'descend');
sortData = cat(dim,pSort,nSort);
end