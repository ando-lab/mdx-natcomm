classdef GridDesigner < util.propertyValueConstructor
    %GRIDDESIGNER
    
    properties
        dgrid = 0.3
        SpaceGroup = symm.SpaceGroup.empty()
        Basis = latt.Basis(1,1,1,90,90,90)
        maxPrimeFactor = 7
        supercell = [1,1,1]
    end
    
    methods
        function obj = GridDesigner(varargin)
            %GRIDDESIGNER
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [LG] = run(obj)
            
            if isempty(obj.SpaceGroup)
                m = [1;1;1];
            else
                m = obj.find_grid_factors();
            end
            
            nx = next_good_dimension(obj,obj.Basis.a/obj.dgrid,m(1));
            ny = next_good_dimension(obj,obj.Basis.b/obj.dgrid,m(2));
            nz = next_good_dimension(obj,obj.Basis.c/obj.dgrid,m(3));
            
            % make the lattice grid object
            P = latt.PeriodicGrid([nx,ny,nz].*obj.supercell,[0,0,0],obj.supercell);
            LG = latt.LatticeGrid(P,obj.Basis);
            
        end
        
        
        
        function m = find_grid_factors(obj)
            
            ops = obj.SpaceGroup.generalPositions;
            
            r = [1,1,1]';
            
            for j = 1:numel(ops)
                t = ops(j).t;
                isLess = t~=0 & t < r;
                r(isLess) = t(isLess);
            end
            
            m = round(1./r); % round just in case
            
        end
        
        
        function n = next_good_dimension(obj,nMin,requiredFactor)
            
            if isempty(obj.maxPrimeFactor)
                f = Inf;
            else
                assert(numel(factor(obj.maxPrimeFactor))==1,'maxPrimeFactor must be prime!');
                f = obj.maxPrimeFactor;
            end
            
            n = floor(nMin);
            while mod(n,requiredFactor) ~= 0 || max(factor(n)) > f
                n = n + 1;
            end
            
        end
        
        
        
        
    end
    
    methods(Static)
        function s = slimit(LG)
            B = LG.Basis.orthogonalizationMatrix;
            P = LG.PeriodicGrid;
            [h,k,l] = P.invert.ind2frac([1,P.N(1)],[1,P.N(2)],[1,P.N(3)]);
            hmax = min(abs(h));
            kmax = min(abs(k));
            lmax = min(abs(l));
            
            [s] = h2slimit(B,hmax,kmax,lmax);
        end
        
        function [Gx,Gy,Gz] = orthoslice(G,px,py,pz)
            
            if isa(G,'latt.PeriodicGrid')
                [n1,n2,n3] = G.frac2ind(px,py,pz);
                
                Gx = G; Gx.P(1) = Gx.P(1)/Gx.N(1); Gx.N(1) = 1;
                Gy = G; Gy.P(2) = Gy.P(2)/Gy.N(2); Gy.N(2) = 1;
                Gz = G; Gz.P(3) = Gz.P(3)/Gz.N(3); Gz.N(3) = 1;
                
                [Gx.ori(1),Gy.ori(2),Gz.ori(3)] = G.ind2frac(n1,n2,n3);
                
            elseif isa(G,'latt.LatticeGrid')
                [f1,f2,f3] = G.Basis.lab2frac(px,py,pz);
                P = G.PeriodicGrid;
                [Px,Py,Pz] = proc.script.GridDesigner.orthoslice(P,f1,f2,f3);
                Gx = latt.LatticeGrid(Px,G.Basis);
                Gy = latt.LatticeGrid(Py,G.Basis);
                Gz = latt.LatticeGrid(Pz,G.Basis);
            else
                error('G should be latt.PeriodicGrid or latt.LatticeGrid');
            end
            
        end
        
    end
    
end



% function [hmax,kmax,lmax] = s2hlimit(B,smax)
% % B is the orthogonalization matrix of the unit cell
% R = inv(B)'; % <-- ortho matrix of the reciprocal unit cell
% 
% v3 = cross(R(:,1),R(:,2));
% v3 = v3/norm(v3);
% 
% v1 = cross(R(:,2),R(:,3));
% v1 = v1/norm(v1);
% 
% v2 = cross(R(:,3),R(:,1));
% v2 = v2/norm(v2);
% 
% hr = [1,0,0]*B'*v1;
% kr = [0,1,0]*B'*v2;
% lr = [0,0,1]*B'*v3;
% 
% %smax = min([hmax/hr,kmax/kr,lmax/lr]);
% hmax = hr*smax;
% kmax = kr*smax;
% lmax = lr*smax;
% 
% end

function [s] = h2slimit(B,hmax,kmax,lmax)
% B is the orthogonalization matrix of the unit cell
R = inv(B)'; % <-- ortho matrix of the reciprocal unit cell

v3 = cross(R(:,1),R(:,2));
v3 = v3/norm(v3);

v1 = cross(R(:,2),R(:,3));
v1 = v1/norm(v1);

v2 = cross(R(:,3),R(:,1));
v2 = v2/norm(v2);

hr = [1,0,0]*B'*v1;
kr = [0,1,0]*B'*v2;
lr = [0,0,1]*B'*v3;

s = min([hmax/hr,kmax/kr,lmax/lr]);

end
