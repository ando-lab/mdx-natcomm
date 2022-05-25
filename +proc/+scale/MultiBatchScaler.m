classdef MultiBatchScaler < handle
    % fit scaling models to multiple batches
    properties
        I = []
        sigma = []
        ih = [] % index
        ihmax = [] % maximum value of ih
        Model = proc.scale.ScaleFactors.empty();
    end
    methods
        function obj = MultiBatchScaler(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end
        function [scale,offset] = total(obj)
            % Imeas = a*d*(b*I + c)
            %       = scale*(I + offset)
            %   where scale = a*d*b, offset = c/b

            aa = obj.Model.aVal;
            bb = obj.Model.bVal;
            cc = obj.Model.cVal;
            dd = obj.Model.dVal;
            scale = aa.*bb.*dd;
            offset = cc./bb;
        end
        function [Isc,sigmasc] = applyModel(obj)
            % Model:
            %   I = a*d*(b*Isc + c)
            % Solve for Isc:
            %   Isc = I/(a*d*b) - c/b = I/scale - offset
            %   where scale = a*d*b, offset = c/b

            if isempty(obj.sigma)
                thisSigma = 1;
            else
                thisSigma = obj.sigma;
            end
            [scale,offset] = obj.total();
            Isc = obj.I./scale - offset;
            sigmasc = thisSigma./scale;
        end

        function [Ipred] = predict(obj,Imerge)
            % Imeas = a*d*(b*I + c)
            %       = scale*(I + offset)
            [scale,offset] = obj.total();
            Ipred = scale.*(Imerge(obj.ih) + offset);
        end

        function [Im,sigmam] = merge(obj)
            n = obj(1).ihmax;
            wm = zeros(n,1);
            Im = zeros(n,1);
            for j=1:length(obj)
                [Isc,sigmasc] = obj(j).applyModel();
                w = 1./sigmasc.^2;
                wm = wm + accumarray(obj(j).ih,w,[n,1]);
                Im = Im + accumarray(obj(j).ih,Isc.*w,[n,1]);
            end
            wm(wm==0) = Inf; # fix divide by zero bug when all bins are empty
            Im = Im./wm;
            sigmam = 1./sqrt(wm);
        end

        function x2 = bFit(obj,Imerge,lambda)
            if length(obj)>1
                x2 = zeros(1,length(obj));
                for j=1:length(obj)
                    x2(j) = bFit(obj(j),Imerge,lambda);
                end
            else
                [AA,BB,Ay,yy] = obj.bFitProblem(Imerge);
                b = (AA + lambda*(BB))\Ay;
                obj.Model.b = b;

                x2 = b'*AA*b - 2*b'*Ay + yy;
                x2 = x2/numel(obj.I); % reduced chi squared
            end
        end

        function [AA,BB,Ay,yy] = bFitProblem(obj,Imerge)
            [A0,B] = obj.Model.bCalc();
            b0 = obj.Model.b; % save current model
            try
                obj.Model.b = []; % sets to default of 1

                [Isc,sigmasc] = obj.applyModel();
                A = diag(sparse(Imerge(obj.ih)./sigmasc))*A0;
                y = Isc./sigmasc;

                AA = A'*A;
                BB = B'*B;
                Ay = A'*y;
                yy = sum(y.^2);
                obj.Model.b = b0; % put back
            catch EM % just in case...
                obj.Model.b = b0; % put back
                rethrow(EM);
            end
        end

        function [diagAA,diagBB] = bCalcDiag(obj,Imerge)
            % the matrix diagonals can be used to estimate a good value for
            % the regularization parameter.
            diagAA = [];
            diagBB = [];
            if length(obj)>1
                diagAA = [];
                diagBB = [];
                for j=1:length(obj)
                    [thisAA,thisBB] = obj(j).bCalcDiag(Imerge);
                    diagAA = [diagAA; thisAA];
                    diagBB = [diagBB; thisBB];
                end
            else % length(obj) == 1
                [AA,BB] = obj.bFitProblem(Imerge);
                diagAA = full(diag(AA));
                diagBB = full(diag(BB));
            end
        end

        function x2 = dFit(obj,Imerge,lambda)
            % fit globally
            n = prod(obj(1).Model.szd);
            AA = sparse(n,n);
            Ay = sparse(n,1);
            H = sparse(n,n);
            One = sparse(n,1);
            yy = 0;
            nBatches = length(obj);
            for j=1:nBatches
                [thisAA,thisH,thisOne,thisAy,thisyy] = obj(j).dFitProblem(Imerge);
                AA = AA + thisAA;
                Ay = Ay + thisAy;
                H = H + thisH;
                One = One + thisOne;
                yy = yy + thisyy;
            end

            d = (AA + lambda*H)\(Ay + lambda*One);
            N = 0;
            for j=1:nBatches
                obj(j).Model.d = d;
                N = N + numel(obj(j).I);
            end
            x2 = d'*AA*d - 2*d'*Ay + yy;
            x2 = x2/N; % reduced chi squared
        end

        function [diagAA,diagH] = dCalcDiag(obj,Imerge)
            % fit globally
            n = prod(obj(1).Model.szd);
            AA = sparse(n,n);
            H = sparse(n,n);
            nBatches = length(obj);
            for j=1:nBatches
                [thisAA,thisH] = obj(j).dFitProblem(Imerge);
                AA = AA + thisAA;
                H = H + thisH;
            end
            diagAA = full(diag(AA));
            diagH = full(diag(H));
        end

        function [AA,H,One,Ay,yy] = dFitProblem(obj,Imerge)
            A0 = obj.Model.dCalc();
            n = prod(obj.Model.szd);
            H = speye(n);
            One = ones(n,1);
            d0 = obj.Model.d; % save current model
            try
                obj.Model.d = []; % sets to default of 1
                Ipred = obj.predict(Imerge);
                A = diag(sparse(Ipred./obj.sigma))*A0;
                AA = A'*A;
                y = obj.I./obj.sigma;
                yy = sum(y.^2);
                Ay = A'*(y);
                obj.Model.d = d0; % put back
            catch EM % just in case...
                obj.Model.d = d0; % put back
                rethrow(EM);
            end
        end

        function x2 = aFit(obj,Imerge,xyLambda,zLambda)
            if length(obj)>1
                x2 = ones(1,length(obj));
                for j=1:length(obj)
                    x2(j) = aFit(obj(j),Imerge,xyLambda,zLambda);
                end
            else
                [AA,BBxy,BBz,Ay,yy,sz] = aFitProblem(obj,Imerge);
                a = (AA + xyLambda*BBxy + zLambda*BBz)\Ay;
                obj.Model.a = reshape(a,sz);

                x2 = a'*AA*a - 2*a'*Ay + yy;
                x2 = x2/numel(obj.I); % reduced chi squared
            end
        end

        function [AA,BBxy,BBz,Ay,yy,sz] = aFitProblem(obj,Imerge)
            [A0,Bxy,Bz] = obj.Model.aCalc();
            a0 = obj.Model.a; % save current model
            sz = obj.Model.sza;
            try
                obj.Model.a = []; % defaults to a = 1

                Ipred = obj.predict(Imerge);
                y = obj.I./obj.sigma;
                A = diag(sparse(Ipred./obj.sigma))*A0;
                AA = A'*A;
                BBxy = Bxy'*Bxy;
                BBz = Bz'*Bz;
                Ay = A'*y;
                yy = sum(y.^2);
                obj.Model.a = a0; % put back
            catch EM % just in case...
                obj.Model.a = a0; % put back
                rethrow(EM);
            end
        end

        function [diagAA,diagBBxy,diagBBz] = aCalcDiag(obj,Imerge)
            % the matrix diagonals can be used to estimate a good value for
            % the regularization parameter.

            diagAA = [];
            diagBBxy = [];
            diagBBz = [];
            if length(obj)>1
                diagAA = [];
                diagBBxy = [];
                diagBBz = [];
                for j=1:length(obj)
                    [thisAA,thisBBxy,thisBBz] = obj(j).aCalcDiag(Imerge);
                    diagAA = [diagAA; thisAA];
                    diagBBxy = [diagBBxy; thisBBxy];
                    diagBBz = [diagBBz; thisBBz];
                end
            else % length(obj) == 1
                [AA,BBxy,BBz] = obj.aFitProblem(Imerge);
                diagAA = full(diag(AA));
                diagBBxy = full(diag(BBxy));
                diagBBz = full(diag(BBz));
            end
        end

        function x2 = cFit(obj,Imerge,xLambda,yLambda,cLambda)
            if length(obj)>1
                x2 = zeros(1,length(obj));
                for j=1:length(obj)
                    x2(j) = cFit(obj(j),Imerge,xLambda,yLambda,cLambda);
                end
            else
                [AA,BBx,BBy,H,Ay,yy,sz] = obj.cFitProblem(Imerge);
                c = (AA + xLambda*BBx + yLambda*BBy + cLambda*H)\Ay;
                obj.Model.c = reshape(c,sz);

                x2 = c'*AA*c - 2*c'*Ay + yy;
                x2 = x2/numel(obj.I); % reduced chi squared
            end
        end

        function [AA,BBx,BBy,H,Ay,yy,sz] = cFitProblem(obj,Imerge)
            [A,Bx,By,AA] = obj.Model.cCalc();
            sz = obj.Model.szc;
            c0 = obj.Model.c; % save current model
            try
                obj.Model.c = []; % defaults to c = 0
                Ipred = obj.predict(Imerge);
                BBx = Bx'*Bx;
                BBy = By'*By;
                H = speye(prod(sz));
                y = (obj.I-Ipred)./obj.sigma;
                yy = sum(y.^2);
                Ay = A'*y;
                obj.Model.c = c0; % put back
            catch EM % just in case...
                obj.Model.c = c0; % put back
                rethrow(EM);
            end
        end

        function [diagAA,diagBBx,diagBBy,diagH] = cCalcDiag(obj,Imerge)
            % the matrix diagonals can be used to estimate a good value for
            % the regularization parameter.

            if length(obj)>1
                diagAA = [];
                diagBBx = [];
                diagBBy = [];
                diagH = [];
                for j=1:length(obj)
                    [thisAA,thisBBx,thisBBy,thisH] = obj(j).cCalcDiag(Imerge);
                    diagAA = [diagAA; thisAA];
                    diagBBx = [diagBBx; thisBBx];
                    diagBBy = [diagBBy; thisBBy];
                    diagH = [diagH; thisH];
                end
            else % length(obj) == 1
                [AA,BBx,BBy,H] = obj.cFitProblem(Imerge);
                diagAA = full(diag(AA));
                diagBBx = full(diag(BBx));
                diagBBy = full(diag(BBy));
                diagH = full(diag(H));
            end
        end

    end
end
