classdef TinyDetector < geom.Detector
    % 201-by-201 1 mm pixels with no gaps and no chips
    methods
        function obj = TinyDetector(varargin)
            obj@geom.Detector(varargin{:});
            a = defaultArgs();
            n = fieldnames(a);
            for j=1:length(n)
                nj = n{j};
                aj = a.(nj);
                if isempty(obj.(nj))
                    obj.(nj) = aj;
                else
                    if  ~isa(obj.(nj),class(aj)) || ...
                            (ischar(aj) && ~strcmpi(aj,obj.(nj))) ||...
                            (isnumeric(aj) && obj.(nj)~=aj)
                        warning(['argument ''%s'' incompatible with',...
                            'TinyDetector and may cause problems'],nj);
                    end
                end
            end
        end
    end
    methods (Static=true)
        
        function ind = chipIndex(ix,iy,mode)
            % This mimics the behavior of geom.Pilatus6M, but always
            % returns an index of 1
            if nargin==0
                mode = 'all';
            elseif nargin==2
                mode = 'pixels';
            end
            switch lower(mode)
                case 'all'
                    ind = ones(201,201);                    
                case 'pixels'
                    ind = ones(size(ix));
                case 'nearest'
                    ind = ones(size(ix));
                otherwise
                    error('did not recognize mode %s',mode);
            end
        end
    end
end

function args = defaultArgs()
args = struct('name','TinyDetector','nx',201,'ny',201,'qx',1,'qy',1);
end