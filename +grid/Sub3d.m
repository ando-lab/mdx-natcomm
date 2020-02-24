classdef Sub3d < util.propertyValueConstructor
    % Sub3d - a class to calculate reduced indices corresponding to a point's
    % location in space.
    properties
        nbits = 10 % <-- should check that this doesn't overflow
        ndiv % = [1,1,1]
    end
    properties(Dependent = true)
        offset    % = 2^(obj.nbits-1);
        mult      % = 2^obj.nbits;
        arraysize % = [obj.mult^3-1, prod(obj.ndiv)];
    end

    methods
        function obj = Sub3d(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end

        function [wholeIndex,fractionIndex] = getRegion(obj,d,wholeIndex,fractionIndex)
            % d is the number of fractional steps to take in each direction
            % (it must be a whole number)
            %
            % (h-h0)^2/ndiv1^2 + (k-k0)^2/ndiv2^2 + (l-l0)^2/ndiv3^2 <= d^2
            %
            % if only wholeIndex is given, then d is taken to refer to the
            % number of whole steps, rather than the number of fractional
            % steps, returning whole index corresponding to h,k,l with
            % (h-h0)^2 + (k-k0)^2 + (l-l0)^2 <= d^2

            assert(nargin==3 | nargin==4,...
                'incorrect number of input arguments');
            assert(isempty(d) | length(d)==1 | length(d)==3,...
                'd is not of the expected dimensions');

            if nargin==3
                assert(nargout <= 1,'only one output will be assigned');
                [h0,k0,l0] = index2hkl(obj,wholeIndex(:));
                hStep = 1; %take whole steps
                kStep = 1;
                lStep = 1;
            else
                [h0,k0,l0] = index2hkl(obj,wholeIndex(:),fractionIndex(:));
                hStep = 1/obj.ndiv(1); %take fractional steps
                kStep = 1/obj.ndiv(2);
                lStep = 1/obj.ndiv(3);
            end

            if isempty(d)
                d = 1; % default = cartesian nearest neighbors
                dh = [0,-1, 1, 0, 0, 0, 0]*hStep;
                dk = [0, 0, 0,-1, 1, 0, 0]*kStep;
                dl = [0, 0, 0, 0, 0,-1, 1]*lStep;
            elseif length(d) == 3
                % find neighbors in a cartesian grid
                [dh,dk,dl] = ndgrid(-d(1):d(1),-d(2):d(2),-d(3):d(3));
                dh = dh(:)'*hStep;
                dk = dk(:)'*kStep;
                dl = dl(:)'*lStep;
            else% length(d)==1
                % find neighbors (ball)
                [dh,dk,dl] = ndgrid(-d:d,-d:d,-d:d);
                isInside = dh.*dh + dk.*dk + dl.*dl <= d.*d;
                dh = dh(isInside)'*hStep;
                dk = dk(isInside)'*kStep;
                dl = dl(isInside)'*lStep;
            end

            if nargout==2
                [wholeIndex,fractionIndex] = obj.hkl2index(...
                    repmat(h0,1,size(dh,2))+repmat(dh,size(h0,1),1),...
                    repmat(k0,1,size(dk,2))+repmat(dk,size(k0,1),1),...
                    repmat(l0,1,size(dl,2))+repmat(dl,size(l0,1),1));
            else
                
                [wholeIndex] = obj.hkl2index(...
                    repmat(h0,1,size(dh,2))+repmat(dh,size(h0,1),1),...
                    repmat(k0,1,size(dk,2))+repmat(dk,size(k0,1),1),...
                    repmat(l0,1,size(dl,2))+repmat(dl,size(l0,1),1));
            end

        end

        function [ind1,ind2] = regrid(obj,oldGrid,oldInd1,oldInd2)
            [h,k,l] = oldGrid.index2hkl(oldInd1,oldInd2);
            [ind1,ind2] = obj.hkl2index(h,k,l);
        end

        function im = mapArray(obj,ind,val,varargin)
            if nargin==6
                h = varargin{1};
                k = varargin{2};
                l = varargin{3};
                [wholeIndex,fractionIndex] = obj.hkl2index(h,k,l);
            elseif nargin==5
                wholeIndex = varargin{1};
                fractionIndex = varargin{2};
            elseif nargin==4
                [wholeIndex,fractionIndex] = sub2ind(obj.arraysize,varargin{1});
            end

            [mask,rowIndex] = ismember(wholeIndex,ind);
            rowIndex(~mask) = NaN;
            fractionIndex(~mask) = NaN;
            sz = [size(ind,1),obj.arraysize(2)];
            ind1 = sub2ind(sz,rowIndex,fractionIndex);
            isincl = ~isnan(ind1);

            if islogical(val)
                im = false(size(ind1));
            else
                im = NaN*zeros(size(ind1));
            end
            im(isincl) = val(ind1(isincl));
        end
        
        function [h,k,l,dh,dk,dl] = hkl2hkl(obj,varargin)
            % [h,k,l] = hkl2hkl(h,k,l,dh,dk,dl)
            % [h,k,l,dh,dk,dl] = hkl2hkl(h,k,l)
            
            nx = obj.ndiv(1);
            ny = obj.ndiv(2);
            nz = obj.ndiv(3);
            
            hh = varargin{1};
            kk = varargin{2};
            ll = varargin{3};
            
            if length(varargin)==6
                % map d,h,l,dh,dk,dl to h,k,l (fractional)
                assert(nargout <= 3,'[h,k,l] output expected');
                
                dhh = varargin{4};
                dkk = varargin{5};
                dll = varargin{6};

                h = hh + dhh/nx;
                k = kk + dkk/ny;
                l = ll + dll/nz;
                
            elseif length(varargin) == 3
                % map d,h,l,dh,dk,dl to h,k,l (fractional)
                
                h = ceil(hh-0.5);
                k = ceil(kk-0.5);
                l = ceil(ll-0.5);

                dh = ceil(nx*(hh-h+.5)) - 0.5*(nx+1);
                dk = ceil(ny*(kk-k+.5)) - 0.5*(ny+1);
                dl = ceil(nz*(ll-l+.5)) - 0.5*(nz+1);
                
            else
                error('incorrect number of input arguments');
            end
        end

        function [hdec,kdec,ldec,dh,dk,dl] = index2hkl(obj,wholeIndex,fractionIndex)
            m = obj.mult;
            minv = 1/m;
            wholeIndex = wholeIndex(:) - 1; % switch to zero-indexing first!
            ldec = floor(wholeIndex*minv*minv);
            wholeIndex = wholeIndex - ldec*m*m;
            kdec = floor(wholeIndex*minv);
            wholeIndex = wholeIndex - kdec*m;
            hdec = wholeIndex;
            hdec = hdec +1 - obj.offset;
            kdec = kdec +1 - obj.offset;
            ldec = ldec +1 - obj.offset;

            if nargin==3 %&& ~isempty(fractionIndex)
                % convert indices to subscripts of the fractional grid
                n1 = obj.ndiv(1);
                n2 = obj.ndiv(2);
                n3 = obj.ndiv(3);

                ndx = fractionIndex(:) - 1; % switch to zero indexing first!!!
                vi = rem(ndx,n1*n2);
                ix3 = (ndx-vi)/(n1*n2);
                ndx = vi;
                vi = rem(ndx,n1);
                ix2 = (ndx - vi)/n1;
                ix1 = vi;

                if nargout ~= 6
                    hdec = hdec + (ix1 + 0.5)/n1 - 0.5;
                    kdec = kdec + (ix2 + 0.5)/n2 - 0.5;
                    ldec = ldec + (ix3 + 0.5)/n3 - 0.5;
                else % return (integer) fractional coordinates
                    % and keep hdec, kdec, ldec as integers also
                    %
                    % note - only integers if n1, n2, n3 are odd
                    dh = ix1 + 0.5 - 0.5*n1;
                    dk = ix2 + 0.5 - 0.5*n2;
                    dl = ix3 + 0.5 - 0.5*n3;
                end
            end
        end

        function [wholeIndex,fractionIndex] = hkl2index(obj,hdec,kdec,ldec,dh,dk,dl)
            % map the points h,k,l to index of nearest integer value
            % (wholeIndex), and map fractional part to fractional grid index
            % (fractionIndex)
            o = obj.offset;
            m = obj.mult;

            if nargin > 4
                n1 = obj.ndiv(1);
                n2 = obj.ndiv(2);
                n3 = obj.ndiv(3);
                hdec = hdec + dh/n1;
                kdec = kdec + dk/n2;
                ldec = ldec + dl/n3;
            end

            h = ceil(hdec-0.5);
            k = ceil(kdec-0.5);
            l = ceil(ldec-0.5);

            wholeIndex = (h+o) + m*(k+o-1) + m*m*(l+o-1);

            if nargout ==2

                n1 = obj.ndiv(1);
                n2 = obj.ndiv(2);
                n3 = obj.ndiv(3);

                ix1 = ceil(n1*(hdec-h+0.5));
                ix2 = ceil(n2*(kdec-k+0.5));
                ix3 = ceil(n3*(ldec-l+0.5));

                fractionIndex = ix1 + (ix2-1)*n1 + (ix3-1)*n1*n2;
            end
        end

        function val = get.offset(obj)
            val = 2^(obj.nbits-1);
        end
        function val = get.mult(obj)
            val = 2^obj.nbits;
        end
        function val = get.arraysize(obj)
            val = [obj.mult^3-1,prod(obj.ndiv)];
        end
    end
end
