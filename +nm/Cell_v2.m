classdef Cell_v2 < util.propertyValueConstructor
    %CELL - model of unit cell contents
    
    properties
        AsymmetricUnit = nm.Group.empty(); % array of Group objects
        Basis = latt.Basis.empty(); % basis vectors of the unit cell
        UnitCellOperators = symm.SymmetryOperator();
    end
    properties (Dependent = true)
        B
    end
    
    methods
        
        function obj = Cell_v2(varargin)
            %CELL
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function Bval = get.B(obj)
            if isa(obj.Basis,'latt.OrientedBasis')
                Om = obj.Basis.orientationMatrix*obj.Basis.orthogonalizationMatrix;
            else
                Om = obj.Basis.orthogonalizationMatrix;
            end
            Bval = symm.AffineTransformation(Om,[0;0;0]);
            
        end
        
        function P = tl2uxyz(obj)
            
            Ops = obj.UnitCellOperators;
            B = obj.B;
            Binv = inv(B);
            numAtoms = numel([obj.AsymmetricUnit.x]);
            numOps = numel(Ops);
            P0 = tl2uxyz(obj.AsymmetricUnit);
            Pn = cell(1,numOps);
            
            for n = 1:numel(Ops)
                Op = B*Ops(n)*Binv;
                R = sparse(Op.r);
                shiftVec = sparse(n,1,1,numOps,1);
                shiftMat = kron(shiftVec,speye(numAtoms));
                Pn{1,n} = kron(shiftMat,R)*P0;
            end
            P = cell2mat(Pn);
        end
        
        function r = unitCellCoordinates(obj)
            Ops = obj.UnitCellOperators;
            B = obj.B;
            numOps = numel(Ops);
            xn = inv(B)*[obj.AsymmetricUnit.r];
            r = cell(1,numel(Ops));
            for n=1:numOps
                r{n} = B*Ops(n)*xn;
            end
        end
        
    end
end

