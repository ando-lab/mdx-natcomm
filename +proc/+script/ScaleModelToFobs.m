classdef ScaleModelToFobs < util.propertyValueConstructor
    %SCALEMODELTOFOBS
    %
    %
    % Parameterization for aniso U based on crystal system:
    %
    %     latticeSystem      U11    U22    U33    U12    U13    U23
    %     ________________    ___    ___    ___    ___    ___    ___
    %
    %     {'Cubic'       }     1      1      1      0      0      0
    %     {'Hexagonal'   }     2      2      1      2      0      0
    %     {'Monoclinic'  }     1      2      3      0      4      0
    %     {'Orthorhombic'}     1      2      3      0      0      0
    %     {'Rhombohedral'}     2      2      1      2      0      0
    %     {'Tetragonal'  }     1      1      2      0      0      0
    %     {'Triclinic'   }     1      2      3      4      5      6
    %
    %
    properties
        Basis
        SpaceGroup
        ktot = 1; % initial guess
        Bsol = 50; % initial guess
        ksol = -.334; % initial guess
        maxFunEval = 5000;
        FunctionTolerance = 1E-8;
    end
    
    methods
        function obj = ScaleModelToFobs(varargin)
            %ScaleModelToFobs
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function params = run(obj,hklTable)
            % GO!
            
            [hklTable.sx,hklTable.sy,hklTable.sz] = obj.Basis.invert.frac2lab(hklTable.h,hklTable.k,hklTable.l);
            hklTable = hklTable(~isnan(hklTable.Fobs),:); % get rid of NaNs
            
            [numParams,param2struct] = obj.parameterize_model();
            
            initialParams = zeros(1,numParams);
            initialParams(1:3) = [obj.ktot,obj.ksol,obj.Bsol];
            
            Tsol = @(s) latt.Blob(1,0).addB(s.Bsol).rescale(s.ksol);
            Taniso = @(s) latt.Blob(1,0).addU(s.U).rescale(s.ktot);
            
            fmodel = @(s,t) t.Fatom + Tsol(s).scatteringAmplitude(t.sx,t.sy,t.sz).*t.Fsol;
            residfun = @(s,t) (t.Fobs - Taniso(s).scatteringAmplitude(t.sx,t.sy,t.sz).*abs(fmodel(s,t)))./t.sigmaFobs;
            
            opts = optimoptions('lsqnonlin','MaxFunctionEvaluations',obj.maxFunEval,'FunctionTolerance',obj.FunctionTolerance);
            
            solvFit = lsqnonlin(@(v) residfun(param2struct(v/1000),hklTable),1000*initialParams,[],[],opts);
            
            params = param2struct(solvFit/1000);
        end
        
        function [numParams,param2struct] = parameterize_model(obj)
            switch lower(obj.SpaceGroup.Info.latticeSystem)
                case 'cubic'
                    param2struct = @(v) struct(...
                        'ktot',v(1),...
                        'U', [v(4),0,0;0,v(4),0;0,0,v(4)],...
                        'ksol',v(2),...
                        'Bsol',v(3));
                    numParams = 4;
                case {'hexagonal','rhombohedral'}
                    param2struct = @(v) struct(...
                        'ktot',v(1),...
                        'U', [v(4),v(4),0;v(4),v(4),0;0,0,v(5)],...
                        'ksol',v(2),...
                        'Bsol',v(3));
                    numParams = 5;
                case 'monoclinic'
                    param2struct = @(v) struct(...
                        'ktot',v(1),...
                        'U', [v(4),0,v(7);0,v(5),0;v(7),0,v(6)],...
                        'ksol',v(2),...
                        'Bsol',v(3));
                    numParams = 7;
                case 'orthorhombic'
                    param2struct = @(v) struct(...
                        'ktot',v(1),...
                        'U', [v(4),0,0;0,v(5),0;0,0,v(6)],...
                        'ksol',v(2),...
                        'Bsol',v(3));
                    numParams = 6;
                case 'tetragonal'
                    param2struct = @(v) struct(...
                        'ktot',v(1),...
                        'U', [v(4),0,0;0,v(4),0;0,0,v(5)],...
                        'ksol',v(2),...
                        'Bsol',v(3));
                    numParams = 5;
                case 'triclinic'
                    param2struct = @(v) struct(...
                        'ktot',v(1),...
                        'U', [v(4),v(7),v(8);v(7),v(5),v(9);v(8),v(9),v(6)],...
                        'ksol',v(2),...
                        'Bsol',v(3));
                    numParams = 9;
                otherwise
                    error('oops');
            end
        end
        
    end
end

%% the following computes the parameterization by brute force for each space group
%
% vars = {'spaceGroup','U11','U22','U33','U12','U13','U23'};
% bigT = table([],[],[],[],[],[],[],'VariableNames',vars);
%
% for j=1:230
%
% SpaceGroup = symm.SpaceGroup(j);
% Ops = SpaceGroup.generalPositions;
%
% syms U11 U22 U33 U12 U13 U23
% U = [U11 U12 U13; U12 U22 U23; U13 U23 U33];
%
% eqns = [];
%
% for n=1:numel(Ops)
%     eq = U == Ops(n).r*U*(Ops(n).r');
%     eqns = [eqns; eq(:)];
% end
%
% [A,b] = equationsToMatrix(eqns,[U11 U22 U33 U12 U13 U23]);
%
% [r,c] = find(null(A));
%
% T = array2table([SpaceGroup.number,accumarray(r,c,[6,1])'],'VariableNames',vars);
%
% bigT = [bigT;T];
%
% end
%
% bigT.latticeSystem = arrayfun(@(n) symm.SpaceGroup(n).Info.latticeSystem,bigT.spaceGroup,'Uni',0);
%
% disp(unique(bigT(:,{'pointGroup','U11','U22','U33','U12','U13','U23'}),'rows'))