function [data,colNames] = read_integrate_hkl(fileName)

colNames = {};

% seek end of header
fid = fopen(fileName,'r');

lengthOfHeader = 0;

while ~feof(fid)
    thisLine = fgetl(fid);
    
    isHeaderLine = strncmp(thisLine,'!',1);
    isLastHeaderLine = strncmp(thisLine,'!END_OF_HEADER',14);
    
    if isHeaderLine
        lengthOfHeader = lengthOfHeader + 1;
    else
        error('this should never happen');
    end
    
    if isLastHeaderLine
        break;
    end
    
    if strncmp(thisLine,'!H,K,L',6)
        %  !H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK,CORR,MAXC,
        %  !             XOBS,YOBS,ZOBS,ALF0,BET0,ALF1,BET1,PSI,ISEG
        nextLine = fgetl(fid);
        lengthOfHeader = lengthOfHeader + 1;
        colNames = strtrim(strsplit([thisLine(2:end),nextLine(2:end)],','));
    end
    
end

nCols = length(colNames);

data = fscanf(fid,'%f',[nCols,Inf])';

fclose(fid);

%data = dlmread(fileName,' ',lengthOfHeader,0);

end