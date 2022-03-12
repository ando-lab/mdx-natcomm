classdef LatticeDynamics < util.propertyValueConstructor
    %LATTICEDYNAMICS vibrational dynamics using Born/Von-Karman method
    
    properties
        V % cell array of Hessian matrices for 3x3x3 array of nearest neighbors
        supercell = [1,1,1]
        M % matrix of generalized masses
        tol_eigd = [];
    end
    properties(Dependent)
        G_sup
    end
    
    methods
        function obj = LatticeDynamics(varargin)
            %LATTICEDYNAMICSMODEL 
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function Vfull = Vfull(obj)
            G_n = latt.PeriodicGrid([3,3,3],[0,0,0],[1,1,1]).invert;
            G = obj.G_sup;
            Vfull = expand2array(obj.V,G_n,G);
        end
        
        function V = K2V(obj,K)
            V = fourier_transform(K,obj.G_sup.invert);
        end
        
        function K = V2K(obj,V)
            K = inverse_fourier_transform(V,obj.G_sup);
        end
        
        function K = K(obj)
            if all(obj.supercell==1)
                K = obj.Vfull();
            else
                K = inverse_fourier_transform(obj.Vfull(),obj.G_sup);
            end
        end
        
        function Kinv = Kinv(obj)
            Kinv = pseudo_inverse(obj.K(),obj.G_sup,full(obj.M),obj.tol_eigd);
        end
        
        function covMat = cov_unitcell(obj)
            if all(obj.supercell==1)
                covMat = obj.Kinv();
            else
                covMat = mean(obj.Kinv(),[1,2,3]);
                covMat = real(covMat);
            end
            covMat = shiftdim(covMat,3);
        end
        
        function covMat = cov_supercell(obj)
            covMat = fourier_transform(obj.Kinv(),obj.G_sup.invert);
            covMat = real(covMat);
            covMat = shiftdim(covMat,3);
        end
        
        function G_sup = get.G_sup(obj)
            G_sup = latt.PeriodicGrid(obj.supercell,[0,0,0],[1,1,1]).invert;
        end
        
        
%         function I = onePhononScattering(obj,h,k,l,G)
%             Kinv = obj.Kinv;
%             Kinv = shiftdim(Kinv,3);
%             G = cat(2,G{:}); % flatten
%             G_bz = LD.G_sup.invert;
%             [n1,n2,n3] = G_bz.frac2ind(h(:),k(:),l(:));
%             ind = accumarray([n1,n2,n3],1:numel(h),G_bz.N,@(v) {v});
%             
%             sz = size(h);
%             I = zeros(sz);
%             fprintf(1,'accumulating intensitites\n');
%             for j=1:numel(ind)
%                 fprintf(1,'k = %d of %d\n',j,numel(ind));
%                 if isempty(ind{j}), continue; end
%                 Gk = G(ind{j},:);
%                 I(ind{j}) = I(ind{j}) + real(dot(Gk,Gk*conj(Kinv(:,:,j)),2));
%             end
%         end
        
    end
end



function Vexp = expand2array(V,G,G_sup)
[n1,n2,n3] = G.grid();
[n1,n2,n3] = G_sup.frac2ind(n1,n2,n3);
sz = size(V{1});
subs = cell(numel(V),1);
v = cell(numel(V),1);
for j=1:numel(V)
    [r,c,v{j}] = find(V{j});
    n = numel(r);
    l1 = repmat(n1(j),[n,1]);
    l2 = repmat(n2(j),[n,1]);
    l3 = repmat(n3(j),[n,1]);
    subs{j} = [l1,l2,l3,r,c];
end

subs = cell2mat(subs);
v = cell2mat(v);

Vexp = accumarray(subs,v,[G_sup.N,sz]);

end


function Y = fourier_transform(X,R)
    Y = X;
    Y = Y*prod(R.delta); % normalize
    for n=1:3
        Y = fft(Y,[],n);
        Y = fftshift(Y,n);
    end

    % do phase factors
     G = R.invert;
     [x,y,z] = G.ind2frac(1:G.N(1),1:G.N(2),1:G.N(3));
     ph = exp(-2i*pi*R.ori(1)*x(:));
     pk = shiftdim(exp(-2i*pi*R.ori(2)*y(:)),-1);
     pl = shiftdim(exp(-2i*pi*R.ori(3)*z(:)),-2);
     Y = Y.*ph.*pk.*pl; % <-- uses broadcasting
end

function X = inverse_fourier_transform(Y,G)
    X = Y;
    for n=1:3
        X = ifft(X,[],n);
        X = fftshift(X,n);
    end
    
     % do phase factors
     R = G.invert;
     [x,y,z] = R.ind2frac(1:R.N(1),1:R.N(2),1:R.N(3));
     ph = exp(2i*pi*G.ori(1)*x(:));
     pk = shiftdim(exp(2i*pi*G.ori(2)*y(:)),-1);
     pl = shiftdim(exp(2i*pi*G.ori(3)*z(:)),-2);
     X = X.*ph.*pk.*pl; % <-- uses broadcasting
     X = X/prod(R.delta);
end

function Kinv = pseudo_inverse(K,G_sup,M,tol_eigd)
    if nargin < 4 
        tol_eigd = [];
    end
    
    [o1,o2,o3] = G_sup.frac2ind(0,0,0);
    Kinv = zeros(size(K));
    for n1=1:size(K,1)
        for n2=1:size(K,2)
            for n3=1:size(K,3)
                Kn = squeeze(K(n1,n2,n3,:,:));
                if all([n1,n2,n3]==[o1,o2,o3])
                    try
                    [v,d] = eig(0.5*(Kn+Kn'),0.5*(M+M'),'vec');
                    catch EM
                        disp(size(M));
                        rethrow(EM)
                    end
                    if isempty(tol_eigd)
                        d(1:3) = Inf; % by default, remove the first three eigenvalues
                    else
                        d(abs(d)<tol_eigd) = Inf;
                    end
                    Kinv(n1,n2,n3,:,:) = v*diag(1./d)*v';
                    %Kinv(n1,n2,n3,:,:) = pinv(Kn);
                else
                    Kinv(n1,n2,n3,:,:) = inv(Kn);
                end
            end
        end
    end
end