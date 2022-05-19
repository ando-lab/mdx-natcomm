classdef MapInterp < util.propertyValueConstructor
    %MapInterp
    
    properties
        MT % proc.script.MapTools object representing the map to interpolate
        target_supercell = [1,1,1]
        target_slim = 1
    end
    methods
        function obj = MapInterp(varargin)
            %MapInterp
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function [To,MTo,msk] = target_hkl(obj)
            
            MTi = obj.MT;
            assert(MTi.type == "intensity");
            
            UnitCellGrid = proc.script.GridDesigner(...
                'Basis',MTi.Basis.invert,...
                'SpaceGroup',MTi.SpaceGroup,...
                'dgrid',0.5/obj.target_slim,...
                'supercell',obj.target_supercell).run();
            
            MTo = MTi;
            MTo.Grid = UnitCellGrid.PeriodicGrid.invert;
            
            msk = MTo.spherical_mask(obj.target_slim) & MTo.isASU();
            To = MTo.array2table(msk,msk); % target points to interpolate on

        end
        
        function [To,wbsq] = sgolay_filter_at_target(obj,I,sigma,To,sradius,minfraction,resid)
            if nargin < 4 || isempty(To)
                To = obj.target_hkl();
            end
            if nargin < 5 || isempty(sradius)
                V1 = obj.MT.Basis.volume/prod(obj.target_supercell);
                sradius = (V1*(3/4)/pi)^(1/3);
            end
            if nargin < 6 || isempty(minfraction)
                minfraction = 0.5;
            end
            
            if nargin >= 7 && ~isempty(resid)
                % do robust re-fitting with bisquare weights
                [Ifilt,sigmaFilt,wbsq] = sgolay_at_hkl(obj.MT,I,sigma,To,sradius,minfraction,resid);
            else
                [Ifilt,sigmaFilt] = sgolay_at_hkl(obj.MT,I,sigma,To,sradius,minfraction);
            end
            
            sigmaFilt(isnan(Ifilt)) = Inf;
            Ifilt(isnan(Ifilt)) = 0;
            
            To.(4) = Ifilt;
            To.(5) = sigmaFilt;
            To.Properties.VariableNames(4:5) = {'I','sigma'};
            
        end
        
        function I = interpolate_from_table(obj,To)
            [~,MTo] = obj.target_hkl();
            I = MTo.table2array(To,'symexpand');
            [h,k,l] = MTo.Grid.grid();
            GI = griddedInterpolant(h,k,l,I,'cubic');
            msk = obj.MT.spherical_mask(obj.target_slim) & obj.MT.isASU();
            tmp = obj.MT.array2table(msk,msk);
            tmp.I = GI(tmp.h,tmp.k,tmp.l);
            I = obj.MT.table2array(tmp,'symexpand');
        end
        
    end
    
    
end


function [Ifilt,sigmaIfilt,wbsq] = sgolay_at_hkl(MTi,I,sigmaI,To,sradius,minfraction,resid)

use_bisquare_weights = false;
return_bisquare_weights = false;

if nargin >= 7 && ~isempty(resid)
    use_bisquare_weights = true;
    return_bisquare_weights = nargout >= 3;
end

npts = size(To,1);

Ifilt = NaN*ones(npts,1);
sigmaIfilt = NaN*ones(npts,1);

sigmaI(isnan(I)) = Inf;
I(isnan(I)) = 0;

[sx,sy,sz] = MTi.Basis.frac2lab(To.h,To.k,To.l);

if return_bisquare_weights
    w2map = zeros(size(I));
    n2map = zeros(size(I));
    indmap = reshape(1:numel(I),size(I));
end

for j=1:npts
    if mod(j,round(npts/1000))==0
        fprintf(1,'%d of %d (%.1f percent)\n',j,npts,100*j/npts);
    end
    [MTc,cropfun] = MTi.resize('radius',sradius,[sx(j),sy(j),sz(j)]);
    Ic = cropfun(I);
    sigmaIc = cropfun(sigmaI);
    sigmaIc(isnan(Ic)) = Inf;
    Ic(isnan(Ic)) = 0;
    if mean(~isinf(sigmaIc(:))) < minfraction
        continue;
    end
    [h,k,l] = MTc.Grid.grid();
    w = 1./sigmaIc;
    
    if use_bisquare_weights
        w2 = bisquare_weights(cropfun(resid));
        w = w.*w2;
    end
    
    [Ifilt(j),~,sigmaIfilt(j)] = sgolay_fit(h-To.h(j),k-To.k(j),l-To.l(j),Ic,w);
    
    if return_bisquare_weights
        indc = cropfun(indmap);
        n2map(indc) = n2map(indc) + 1;
        w2map(indc) = w2map(indc) + w2;
    end
end

if return_bisquare_weights
    wbsq = w2map./n2map;
end

end


function w_out = bisquare_weights(resid)

w_out = zeros(size(resid));

isIncl = ~isnan(resid);
resid = resid(isIncl);

% bisquare weights
m = mad(resid(:),1);
w = (1-(resid/(6*m)).^2).^2;
isExcl = abs(resid)>(6*m);
w(isExcl) = 0;

w_out(isIncl) = w;
end

function [v0,resid,sigmav0] = sgolay_fit(x,y,z,v,w)

B = quad3model(x,y,z);
Bw = B.*w(:);
AA = Bw'*Bw;
Ab = Bw'*(v(:).*w(:));

AAinv = pinv(AA);
c = AAinv*Ab;
v0 = c(end);
sigmav0 = sqrt(AAinv(end,end));
resid = reshape(v(:) - B*c,size(v));
end


function B = quad3model(x,y,z)

B = ones(numel(x),10);
B(:,1) = x(:).^2;
B(:,2) = y(:).^2;
B(:,3) = z(:).^2;
B(:,4) = x(:).*y(:);
B(:,5) = x(:).*z(:);
B(:,6) = y(:).*z(:);
B(:,7) = x(:);
B(:,8) = y(:);
B(:,9) = z(:);
% B(10,:) = 1;

end