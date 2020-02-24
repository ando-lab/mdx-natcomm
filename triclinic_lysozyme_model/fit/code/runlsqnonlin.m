function [history,xfit,fitinfo] = runlsqnonlin(objfun,x0,xmin,xmax,varargin)

% for example:
%varargin = {'Algorithm','trust-region-reflective',...
%    'MaxFunctionEvaluations',1000,'UseParallel',false};

history = struct('iteration',{},'x',{},'resnorm',{});

opts = optimoptions(@lsqnonlin,varargin{:},'OutputFcn',@outfun);

[xfit,~,~,~,fitinfo] = lsqnonlin(objfun,x0,xmin,xmax,opts);

function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             % do nothing
         case 'iter'
         % Concatenate current point and objective function
                history = [history;...
                    struct('iteration',optimValues.iteration,...
                    'x',x,'resnorm',optimValues.resnorm)];
         case 'done'
             % do nothing
         otherwise
     end
end

end