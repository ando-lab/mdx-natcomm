function [hklMerge,ih] = merge(hklTable,ref,intCol,errCol,varargin)
% MERGE - merge redundant observations by weighted least squares
%
% ref is either a table with columns that will be used to index hklTable,
% or it is a list of column names in hklTable that will be used to generate
% a consensus list using unique:

if nargin < 3 || isempty(intCol)
    intCol = 'I';
end
if nargin < 4 || isempty(errCol)
    errCol = 'sigma';
end
[intCol,intColOut] = parse_column_argument(intCol);
[errCol,errColOut] = parse_column_argument(errCol);

if istable(ref)
    refCols = ref.Properties.VariableNames;
    [~,ih] = ismember(hklTable(:,refCols),ref,'rows');
    hklMerge = ref;
    hMax = size(hklMerge,1);
elseif isnumeric(ref)
    ih = ref;
    hMax = max(ih);
    hklMerge = table;
else % ref is a cell string or an index array
    refCols = ref;
    [hklMerge,~,ih] = unique(hklTable(:,refCols),'rows');
    hMax = size(hklMerge,1);
end


w0 = 1./hklTable.(errCol).^2;
w = accumarray(ih,w0,[hMax,1]);
wm = accumarray(ih,hklTable.(intCol).*w0,[hMax,1])./w;

hklMerge.(intColOut) = wm;
hklMerge.(errColOut) = sqrt(1./w);

for j=1:length(varargin)
    [inCol,outCol] = parse_column_argument(varargin{j});
    wm = accumarray(ih,hklTable.(inCol).*w0,[hMax,1])./w;
    hklMerge.(outCol) = wm;
end

end

function [inCol,outCol] = parse_column_argument(arg)
if iscellstr(arg)
    assert(length(arg)==2,'expected 2 variables in argument');
    inCol = arg{1};
    outCol = arg{2};
else
    inCol = arg;
    outCol = arg;
end
end