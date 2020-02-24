classdef Crystal < util.propertyValueConstructor
    properties
        spaceGroupNumber = 1 % P1 by default
        a % unit cell axis length in Angstrom
        b % unit cell axis length in Angstrom
        c % unit cell axis length in Angstrom
        alpha % in degrees
        beta  % in degrees
        gamma % in degrees
        a_axis % 3-vector, a direction
        b_axis % 3-vector, b direction
        c_axis % 3-vector, c direction
    end
    properties(Dependent)
        UnitCell    % geom.lattice.UnitCell object
        SpaceGroup  % geom.lattice.SpaceGroup object
        isValid     % wrapper for geom.lattice.SpaceGroup.checkProperties
        orthogonalizationMatrix % wrapper for geom.lattice.UnitCell.orthogonalizationMatrix;
        rotationMatrix % matrix rotating crystal to spindle coordinate system
        alignmentError % difference between rotated axes and a_axis, b_axis, ...
    end
    methods
        function obj = Crystal(varargin)
            obj@util.propertyValueConstructor(varargin{:});
        end

        function [hasu,kasu,lasu] = hkl2asu(obj,h,k,l)
            % wrapper function for geom.Scattering.hkl2asu
            [hasu,kasu,lasu] = geom.Scattering.hkl2asu(h,k,l,obj);
        end
        
        function [hnew,knew,lnew,ind] = asu2hkl(obj,h,k,l)
            % wrapper function for geom.Scattering.asu2hkl
            [hnew,knew,lnew,ind] = geom.Scattering.asu2hkl(h,k,l,obj);
        end
        
        function [sx,sy,sz] = hkl2s(obj,h,k,l)
            [sx,sy,sz] = geom.Scattering.hkl2s0(h,k,l,obj);
        end

        function value = get.SpaceGroup(obj)
            value = geom.lattice.SpaceGroup(obj.spaceGroupNumber);
        end

        function value = get.UnitCell(obj)
            value = geom.lattice.UnitCell(obj.a,obj.b,obj.c,...
                obj.alpha,obj.beta,obj.gamma);
        end

        function tf = get.isValid(obj)
            tf = obj.SpaceGroup.checkProperties(obj.a,obj.b,obj.c,...
                obj.alpha,obj.beta,obj.gamma);
        end

        function value = get.orthogonalizationMatrix(obj)
            value = obj.UnitCell.orthogonalizationMatrix;
        end

        function value = get.rotationMatrix(obj)
            targetAxes = [obj.a_axis(:),obj.b_axis(:),obj.c_axis(:)]; %axes are columns
            if isempty(targetAxes)
                value = eye(3);
                return;
            end
            crystalAxes = [obj.UnitCell.a1(:),...
                obj.UnitCell.a2(:),...
                obj.UnitCell.a3(:)];
            [u,~,v] = svd(crystalAxes*targetAxes','econ');
            value = v*u';
            assert(det(value)>0,'coordinate system is not right-handed?');
            
        end

        function value = get.alignmentError(obj)
            targetAxes = [obj.a_axis(:),obj.b_axis(:),obj.c_axis(:)]; %axes are columns
            crystalAxes = [obj.UnitCell.a1(:),...
                obj.UnitCell.a2(:),...
                obj.UnitCell.a3(:)];
            value = obj.rotationMatrix*crystalAxes - targetAxes;
        end
    end
end
