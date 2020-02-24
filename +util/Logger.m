classdef Logger < handle
    %LOGGER extends logging capability of diary()
    properties(SetAccess = immutable)
        fileName % set by constructor
    end
    properties(SetAccess = private)
        isPaused = false; % keeps track of whether the timer is paused
    end
    properties(Access = private)
        t
    end
    properties(Dependent = true)
        isRunning
    end
    methods
        function obj = Logger(fileName)
            % first, check to see if there are other active Loggers
            logTimers = timerfind('UserData',class(obj));
            
            if ~isempty(logTimers)
                error('Found %d active %s object(s). To avoid conflict, only one is allowed a time. Delete the others before creating a new one.',length(logTimers),class(obj));
            end
            
            obj.fileName = fileName;
            obj.t = timer('timerfcn',@updateDiary,...
                'stopfcn',@stopDiary,'period',5,...
                'ExecutionMode','fixedRate',...
                'UserData',class(obj));
        end
        function val = get.isRunning(obj)
            switch lower(obj.t.Running)
                case 'on'
                    val = true;
                case 'off'
                    val = false;
            end
        end
        function start(obj)
            diary(obj.fileName);
            start(obj.t);
        end
        function pause(obj)
            if obj.isRunning
                obj.isPaused = true;
                stop(obj.t);
            end
        end
        function resume(obj)
            if obj.isPaused
                start(obj.t);
                obj.isPaused = false;
            end
        end
        function fprintf(obj,varargin)
            % print something to the log file
            obj.pause;
            fopen(fid,obj.fileName,'a');
            fprintf(fid,varargin{:});
            fclose(fid);
            obj.resume;
        end
        function clear(obj)
            if obj.isRunning
                stop(obj);
                delete(obj.fileName);
                start(obj);
            else
                delete(obj.fileName);
            end
        end
        function stop(obj)
            stop(obj.t);
            obj.isPaused = false; % clears pause also 
        end
        function delete(obj)
            stop(obj.t);
            delete(obj.t);
        end
    end
end

function updateDiary(~,~,~)
diary off;
diary on;
end

function stopDiary(~,~,~)
diary off;
end