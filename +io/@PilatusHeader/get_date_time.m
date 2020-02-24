function value = get_date_time(obj)
date_time_pattern = ['(\d{4}-\d{2}-\d{2}T|\d{4}/\D+/\d{2} )',...
    '\d{2}:\d{2}:\d{2}.\d+'];
dtmatch = regexp(obj.rawheader,date_time_pattern,'match');
if ~isempty(dtmatch)
    value = dtmatch{1};
else
    value = '';
end
end