function [wedge,suggestedValues] = read_integrate_lp(fileName)

fid = fopen(fileName,'r');

xparm_in = read_header_xparm(fid);

j = 1;
wedge = struct('dataRange',{},'table',{},'stats','xparm');

while ~feof(fid)
    newline = fgetl(fid);
    
    if iswedgestart(newline)
        tmp_cell= regexp(newline,...
            '^\s+PROCESSING OF IMAGES\s+(\d+) \.{3}\s+(\d+)\s*$','tokens');
        n1 = str2double(tmp_cell{1}{1});
        n2 = str2double(tmp_cell{1}{2});
        wedge(j).dataRange = [n1,n2];
    end
    
    if istablestart(newline)
        nLines = wedge(j).dataRange(2) - wedge(j).dataRange(1) + 1;
        wedge(j).table = read_table(fid,nLines);
    end
    
    if isstatsstart(newline)
        wedge(j).stats = read_stats(fid);
        startingFrame = wedge(j).dataRange(1);
        wedge(j).xparm = update_xparm(wedge(j).stats,startingFrame,xparm_in);
        % remove fields which are redundant with xparm.
        wedge(j).stats = rmfield(wedge(j).stats,...
            {'unit_cell','a_axis','b_axis','c_axis','rotation_axis',...
            'f','direct_beam_coords','detector_origin_pixels'});
    end
    
    if iswedgeend(newline)
        j = j + 1;
    end
    
    if issuggestedvalues(newline)
        newline = fgetl(fid);
        tmp1 = sscanf(newline,' BEAM_DIVERGENCE= %f BEAM_DIVERGENCE_E.S.D.= %f');
        newline = fgetl(fid);
        tmp2 = sscanf(newline,' REFLECTING_RANGE= %f REFLECTING_RANGE_E.S.D.= %f');
        suggestedValues.beam_divergence = tmp1(1);
        suggestedValues.beam_divergence_esd = tmp1(2);
        suggestedValues.reflecting_range = tmp2(1);
        suggestedValues.reflecting_range_esd = tmp2(2);
    end
    
end
fclose(fid);

end


% helper functions

function tf = istablestart(newline)
tf = ~isempty(regexp(newline,'^\s*IMAGE\s+IER\s+SCALE','once'));
end

function tf = isstatsstart(newline)
tf = ~isempty(regexp(newline,...
    '^\s*\d+\s+OUT OF\s+\d+\s+REFLECTIONS ACCEPTED FOR REFINEMENT',...
    'once'));
end

function tf = iswedgestart(newline)
tf = ~isempty(regexp(newline,'^\s+PROCESSING OF IMAGES','once'));
end

function tf = iswedgeend(newline)
tf = ~isempty(regexp(newline,'^\s*\*{5} AVERAGE','once'));
end

function tf = issuggestedvalues(newline)
tf = ~isempty(regexp(newline,'^\s*\*{5} SUGGESTED VALUES','once'));
end

