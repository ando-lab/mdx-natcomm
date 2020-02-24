function [Cryst1,Atom,Hetatm,Anisou,remarkLines] = read(fileName)

Cryst1 = struct('a',{},'b',{},'c',{},'alpha',{},'beta',{},'gamma',{},'sGroup',{},'z',{});
Atom = table();
Hetatm = table();
Anisou = table();

pdbStr = fileread(fileName);

atomExpr = '(?<=(\n|^))ATOM\s\s.*?(?=(\n|$))';
hetatmExpr = '(?<=(\n|^))HETATM.*?(?=(\n|$))';
anisouExpr = '(?<=(\n|^))ANISOU.*?(?=(\n|$))';
cryst1Expr = '(?<=(\n|^))CRYST1.*?(?=(\n|$))';
remarkExpr = '(?<=(\n|%))REMARK.*?(?=(\n|$))';

atomLines = regexp(pdbStr,atomExpr,'match');
hetatmLines = regexp(pdbStr,hetatmExpr,'match');
anisouLines = regexp(pdbStr,anisouExpr,'match');
cryst1Lines = regexp(pdbStr,cryst1Expr,'match');
remarkLines = regexp(pdbStr,remarkExpr,'match');

% parse cryst1 record
if ~isempty(cryst1Lines)
    assert(numel(cryst1Lines)==1,'more than one Cryst1 record?');
    
    % the following ensures at least 80 characters per line
    firstLine = repmat(' ',1,80); % add blank line of 80 characters
    allLines = char([{firstLine},cryst1Lines]);
    allLines = allLines(2:end,:); % remove blank line at begining
        
    Cryst1 = struct('a',str2double(allLines(7:15)),...
        'b',str2double(allLines(16:24)),...
        'c',str2double(allLines(25:33)),...
        'alpha',str2double(allLines(34:40)),...
        'beta',str2double(allLines(41:47)),...
        'gamma',str2double(allLines(48:54)),...
        'sGroup',strtrim(allLines(56:66)),...
        'z',str2double(allLines(67:end)));
end

% parse atom records
if ~isempty(atomLines)
    
    % the following ensures at least 80 characters per line
    firstLine = repmat(' ',1,80); % add blank line of 80 characters
    allLines = char([{firstLine},atomLines]);
    allLines = allLines(2:end,:); % remove blank line at begining
    
    % assign records
    serial = str2num(allLines(:,7:11));
    name = strtrim(cellstr(allLines(:,13:16)));
    altLoc = strtrim(cellstr(allLines(:,17)));
    resName = strtrim(cellstr(allLines(:,18:20)));
    chainID = strtrim(cellstr(allLines(:,22)));
    resSeq = str2num(allLines(:,23:26));
    iCode = strtrim(cellstr(allLines(:,27)));
    x = str2num(allLines(:,31:38));
    y = str2num(allLines(:,39:46));
    z = str2num(allLines(:,47:54));
    occupancy = str2num(allLines(:,55:60));
    tempFactor = str2num(allLines(:,61:66));
    element = strtrim(cellstr(allLines(:,77:78)));
    charge = strtrim(cellstr(allLines(:,79:80)));
    
    Atom = table(serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor,element,charge);
end

% parse hetatm records
if ~isempty(hetatmLines)
    
    % the following ensures at least 80 characters per line
    firstLine = repmat(' ',1,80); % add blank line of 80 characters
    allLines = char([{firstLine},hetatmLines]);
    allLines = allLines(2:end,:); % remove blank line at begining
    
    % assign records
    serial = str2num(allLines(:,7:11));
    name = strtrim(cellstr(allLines(:,13:16)));
    altLoc = strtrim(cellstr(allLines(:,17)));
    resName = strtrim(cellstr(allLines(:,18:20)));
    chainID = strtrim(cellstr(allLines(:,22)));
    resSeq = str2num(allLines(:,23:26));
    iCode = strtrim(cellstr(allLines(:,27)));
    x = str2num(allLines(:,31:38));
    y = str2num(allLines(:,39:46));
    z = str2num(allLines(:,47:54));
    occupancy = str2num(allLines(:,55:60));
    tempFactor = str2num(allLines(:,61:66));
    element = strtrim(cellstr(allLines(:,77:78)));
    charge = strtrim(cellstr(allLines(:,79:80)));
    
    Hetatm = table(serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor,element,charge);
end

% parse anisou lines
if ~isempty(anisouLines)
    
    % the following ensures at least 80 characters per line
    firstLine = repmat(' ',1,80); % add blank line of 80 characters
    allLines = char([{firstLine},anisouLines]);
    allLines = allLines(2:end,:); % remove blank line at begining
    
    % assign records
    serial = str2num(allLines(:,7:11));
    name = strtrim(cellstr(allLines(:,13:16)));
    altLoc = strtrim(cellstr(allLines(:,17)));
    resName = strtrim(cellstr(allLines(:,18:20)));
    chainID = strtrim(cellstr(allLines(:,22)));
    resSeq = str2num(allLines(:,23:26));
    iCode = strtrim(cellstr(allLines(:,27)));
    u00 = str2num(allLines(:,29:35));
    u11 = str2num(allLines(:,36:42));
    u22 = str2num(allLines(:,43:49));
    u01 = str2num(allLines(:,50:56));
    u02 = str2num(allLines(:,57:63));
    u12 = str2num(allLines(:,64:70));
    element = strtrim(cellstr(allLines(:,77:78)));
    charge = strtrim(cellstr(allLines(:,79:80)));
    
    Anisou = table(serial,name,altLoc,resName,chainID,resSeq,iCode,u00,u11,u22,u01,u02,u12,element,charge);

end

end


function v = str2num(A) % override str2num

    v = zeros(size(A,1),1);
    for j=1:size(v,1)
        v(j) = str2double(A(j,:));
    end
end