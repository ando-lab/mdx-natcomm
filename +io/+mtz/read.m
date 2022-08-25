function [ hklArray, colNames, crystal ] = read( fileName )
% IO.MTZ.READ - function to read reflection data from an mtz file 
%
% Usage:
% [ hklArray, colNames, crystal ] = read( fileName )
%
% Returns:
%   hklArray - an array containing the mtz file data
%   colNames - a cell array of strings with the names of each column
%    crystal - a struct object representing the unit cell with fields:
%              spaceGroupNumber, a, b, c, alpha, beta, gamma
%
% Resources:
%   ccp4 mtz spec: http://www.ccp4.ac.uk/html/mtzformat.html

% open the file and read some header bytes
fid = fopen(fileName,'r');

try
    
A1 = fread(fid,[1,4],'char*1=>char'); % now we're 4 bytes in. char(A1) should be 'MTZ '
assert(strcmp(A1,'MTZ '),'expected ''MTZ '' at the start of the file');
A2 = fread(fid,[1,1],'int32'); % now we're 8 bytes in. A2 is start of the header records

% read the text header at the end of the file
fseek(fid,A2*4,'bof'); % go to position A2
header = fread(fid,[1,Inf],'char*1=>char');

% parse header to determine the number of columns
ncol = regexp(header,'NCOL\s+(?<ncols>\d+)\s+(?<nrows>\d+)','names');
ncols = str2double(ncol.ncols);
nrows = str2double(ncol.nrows);

% calculate the number of rows
% nrows = (A2 - 21)/ncols; % the number of rows to read in

% read the byte data into a table
fseek(fid,80,'bof'); % start reading at row 21
hklArray = fread(fid,[ncols,nrows],'real*4=>double')';

% close the file
fclose(fid);

catch EM
    fclose(fid); % just in case it throws an error before I can close the file
    rethrow(EM);
end

colNames = regexp(header,'COLUMN\s+(\S+)','tokens');
colNames = [colNames{:}];
assert(length(colNames)==ncols,'number of column names did not match number of columns');

cellInfo = regexp(header,'CELL ([\d\s\.]+)','match');
cellInfo = cellInfo{1}; % just use the first appearance
cellParams = sscanf(cellInfo,'CELL %f %f %f %f %f %f',[1,6]);

symInfo = regexp(header,'SYMINF\s+(?<numSymOps>\d)+\s+(?<numPrimOps>\d)+\s+(?<latticeType>[A-Z]+)\s+(?<spaceGroupNumber>\d+)\s+\''(?<spaceGroup>.*)\''\s+(?<pointGroup>\S+)','names');

%Crystal = geom.Crystal('spaceGroupNumber',spaceGroupNumber,...
%    'a',cellParams(1),'b',cellParams(2),'c',cellParams(3),...
%    'alpha',cellParams(4),'beta',cellParams(5),'gamma',cellParams(6));

crystal = struct('spaceGroupNumber',str2double(symInfo.spaceGroupNumber),...
    'a',cellParams(1),'b',cellParams(2),'c',cellParams(3),...
    'alpha',cellParams(4),'beta',cellParams(5),'gamma',cellParams(6));

end

