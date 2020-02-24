classdef LatticeModel
    %% LATTICEMODEL - calculate normal modes of a crystalline elastic network
    
    properties
        M % mass matrix
        K % cell array of hessian matrices
        N = [1,1,1]; % number of unit cells
        minimumOmega = 1E-6; % smallest value of omega that is not set to zero
    end
    methods
        function obj = LatticeModel(M,K,N)%,SG)
            if nargin >= 1
                obj.M = M;
            end
            if nargin >= 2
                obj.K = K;
            end
            if nargin >= 3
                obj.N = N;
            end
        end
        
        function Dk = dynamicalMatrix(obj,k1,k2,k3)
            [n1,n2,n3] = ndgrid(-1:1,-1:1,-1:1);
            n123 = [n1(:),n2(:),n3(:)];
            Dk = zeros(size(obj.M));
            for n = 1:size(n123,1)
                rcell = n123(n,:);
                phi = exp(2i*pi*(k1*rcell(1) + k2*rcell(2) + k3*rcell(3)));
                Dk = Dk + (obj.K{n}')*phi;
            end
        end
        
        function modeArray = normalModes(obj)
            L = chol(obj.M,'lower');
            Linv = full(inv(L));
            
            P = latt.PeriodicGrid(obj.N,-floor(obj.N/2)./obj.N,[1,1,1]);
            [n1,n2,n3] = ndgrid(1:obj.N(1),1:obj.N(2),1:obj.N(3));
            [k1,k2,k3] = P.ind2frac(n1,n2,n3);
            
            modeArray = struct('k',{},'v',{},'omega',{},'d',{});
            
            for n=1:prod(obj.N)
                Dk = obj.dynamicalMatrix(k1(n),k2(n),k3(n));
                A = Linv*Dk*Linv';
                A = 0.5*A + 0.5*A'; % force symmetric so v*v' is identity
                [v,d] = eig(A);
                omega = real(sqrt(diag(d)));
                
                iszero = omega<obj.minimumOmega; % 1E-7
                omega(iszero) = Inf;
                
                modeArray(n).k = [k1(n),k2(n),k3(n)];
                modeArray(n).v = v;
                modeArray(n).omega = omega;
                modeArray(n).d =d; % save for DEBUGGING
                modeArray(n).Dk = Dk; % save for DEBUGGING
                modeArray(n).Linv = Linv; % save for DEBUGGING
            end
            modeArray = reshape(modeArray,obj.N);
        end
        
        function covMat = covarianceMatrix(obj,modeArray,neighborList)
            
            numCoords = size(obj.M,1);
            nmodes = numel(modeArray);
            covMat = zeros(numCoords,numCoords,size(neighborList,1));
            
            L = chol(obj.M,'lower');
            Linv = full(inv(L));
            
            for k=1:nmodes
                dinv = diag(1./modeArray(k).omega.^2);
                Kinv = Linv'*modeArray(k).v*dinv*modeArray(k).v'*Linv;
                k123 = modeArray(k).k;
                for c=1:size(neighborList,1)
                    rcell = neighborList(c,:);
                    phi = exp(-2i*pi*(k123(1)*rcell(1) + k123(2)*rcell(2) + k123(3)*rcell(3))); % from careful examination of derivation
                    covMat(:,:,c) = covMat(:,:,c) + Kinv*phi;
                end
            end
            
            covMat = real(covMat)/nmodes; 
        end
        
        function [I,Kinv] = onePhononScattering(obj,modeArray,qdotF1)
            I = cell(size(modeArray));
            
            % 1) calculate inverse of spring constants matrix for each k
            L = chol(obj.M,'lower');
            Linv = inv(L);
            Kinv = cell(size(modeArray));
            for k=1:numel(modeArray)
                v = modeArray(k).v;
                omega = modeArray(k).omega;
                Dinv = v*diag(1./omega.^2)*v';
                Kinv{k} = Linv'*Dinv*Linv;
            end
            
            % 2) calculate intensities
            for k=1:numel(modeArray)
                % note: the dot function takes the complex conjugate
                I{k} = real(dot(qdotF1{k}*Kinv{k},qdotF1{k},2));
            end
            
        end
        
    end
end