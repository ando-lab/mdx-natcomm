function [tf, errorMessage] = rescale(varargin)

options = struct(...
    'workingDirectory','./',...
    'unitCellInventoryIn','unitCellInventory.mat',...
    'mergeIn','merge.mat',...
    'mergeBraggIn','mergeBragg.mat',...
    'mergeOut','mergeAbsolute.mat',... % bragg and diffuse together
    'scaleIn','scale.mat',...
    'scaleOut','scaleAbsolute.mat',...
    'matOut','rescale.mat',... % stats about scaling
    'logOut','rescale.log',...
    'smax',1,...
    'npts',501,...
    'scutoff',0.7);

BatchProcess = proc.Batch('proc.Batch.rescale',options,varargin{:});
options = BatchProcess.options; % post-processed options

try % START MAIN SCRIPT
    
    BatchProcess.start();
    
    % PART I: data preprocessing
    %
    % PART II: scattering prediction
    %
    % PART III: find scale factor
    %
    % PART IV: apply scaling corrections
    %
    % PART V: make new scaling model
    
    [Molecules,occupancies,Crystal] = BatchProcess.readFromMatFile(...
        options.unitCellInventoryIn,'Molecules','occupancies','Crystal');
    
    % PART I: data preprocessing
    
    CT = proc.script.CombineTables(...
        'postfix2','Bragg',...
        'hklcols1',{'hasu','kasu','lasu','I','sigma'},...
        'hklcols2',{'hasu','kasu','lasu','I','sigma'},...
        'map2asu',true,...
        'Crystal',Crystal);
    
    hklMerge = BatchProcess.readFromMatFile(options.mergeIn,'hklMerge');
    hklMergeBragg = BatchProcess.readFromMatFile(options.mergeBraggIn,'hklMerge');
    
    % combine diffuse and Bragg into one table
    [hklTable,ia,ib] = CT.run(hklMerge,hklMergeBragg);
    
    % get rid of rows where diffuse data is missing
    hklTable = hklTable(ia,:); % same order as mergeDiffuse.hklMerge
    
    % set missing Bragg data to zero
    isMissing = isnan(hklTable.IBragg);
    hklTable.IBragg(isMissing) = 0;
    hklTable.sigmaBragg(isMissing) = Inf;
    
    ET2A = proc.script.ExpandTableToArray(...
        'hklcols',{'hasu','kasu','lasu','I','IBragg'},...
        'SpaceGroup',symm.SpaceGroup(Crystal.spaceGroupNumber),...
        'symexpand',true,...
        'ndiv',[1,1,1]);
    
    [P,I,IBragg] = ET2A.run(hklTable);
    
    SVR = proc.script.StatisticsVsRadius(...
        'edges',sqrt(linspace(0,options.smax.^2,options.npts)),... % d3s grows linearly with |s| rather than as |s|^2
        'Basis',latt.Basis(Crystal.a,Crystal.b,Crystal.c,Crystal.alpha,Crystal.beta,Crystal.gamma).invert,...
        'PeriodicGrid',P);
    
    % calculate statistics vs radius
    sr = SVR.run(IBragg,I); % <- left off here
    
    % get some variables
    sEdges = SVR.edges(:);
    sCenters = sr.r;
    d3s = sr.d3r;
    IBraggAv = sr.av1;
    IDiffuseAv = sr.av2;
    nObs = sr.n;
    
    % PART II: scattering prediction
    
    Icoh = 0;
    Iincoh = 0;
    Ibond = 0;
    
    for j=1:numel(Molecules)
        [a,b,c] = Molecules(j).scatteringFactors(sCenters);
        Icoh = Icoh + a*occupancies(j);
        Iincoh = Iincoh + b*occupancies(j);
        Ibond = Ibond + c*occupancies(j);
    end
    
    %clear a b c j
    
    nElectrons = sum([Molecules.nElectrons].*occupancies);
    
    Vc = Crystal.UnitCell.vCell;
    
    % PART III: find scale factor
    
    Ipred = -nElectrons^2 + Vc*cumsum(d3s.*(Icoh + Ibond + Iincoh));
    Ipred0 = -nElectrons^2 + Vc*cumsum(d3s.*(Icoh + Iincoh));
    Iobs = Vc*cumsum(d3s.*(IDiffuseAv + IBraggAv));
    sCutoff = sEdges(2:end);
    
    % find scale factor
    ind = find(sCutoff <= options.scutoff,1,'last');
    
    scaleFactor = Ipred(ind)/Iobs(ind);
    scaleFactor0 = Ipred0(ind)/Iobs(ind);
    
    cumulativeIntensity = table(sCutoff,Iobs,Ipred,Ipred0);
    meanIntensity = table(sCenters,IDiffuseAv,IBraggAv,nObs,d3s,Icoh,Ibond,Iincoh,...
        'VariableNames',{'s','Idiffuse','IBragg','nObs','d3s','Icoh','Ibond','Iincoh'});
    meanIntensityAbsolute = table(sCenters,IDiffuseAv*scaleFactor - Iincoh,IBraggAv*scaleFactor,Iincoh,...
        'VariableNames',{'s','Idiffuse','IBragg','Iincoh'});
    
    % save scaling information in output mat file
    BatchProcess.saveToMatFile(options.matOut,...
        'options',options,...
        'scaleFactor',scaleFactor,...
        'scaleFactor0',scaleFactor0,...
        'nElectrons',nElectrons,...
        'Vc',Vc,...
        'cumulativeIntensity',cumulativeIntensity,...
        'meanIntensity',meanIntensity,...
        'meanIntensityAbsolute',meanIntensityAbsolute);
    
    % PART IV: apply scaling corrections
    
    Iincoh = 0;
    
    [sx,sy,sz] = Crystal.hkl2s(hklTable.hasu,hklTable.kasu,hklTable.lasu);
    hklTable.s = sqrt(sx.^2 + sy.^2 + sz.^2);
    
    for j=1:numel(Molecules)
        [~,b] = Molecules(j).scatteringFactors(hklTable.s,'interp');
        Iincoh = Iincoh + b*occupancies(j);
    end
    
    % apply scales
    hklTable.I = hklTable.I*scaleFactor;
    hklTable.sigma = hklTable.sigma*scaleFactor;
    hklTable.IBragg = hklTable.IBragg*scaleFactor;
    hklTable.sigmaBragg = hklTable.sigmaBragg*scaleFactor;
    
    % subtract incoherent scattering
    hklTable.I = hklTable.I - Iincoh;
    hklTable.Iincoh = Iincoh;
    
    % save scaling information in output mat file
    BatchProcess.saveToMatFile(options.mergeOut,...
        'hklMerge',hklTable,...
        'Crystal',Crystal,...
        'scaleFactor',scaleFactor);
    
    
    % % PART V: make new scaling model
    
    ScalingModel = BatchProcess.readFromMatFile(options.scaleIn,'ScalingModel');
    
    
    % Scaling model:
    %   Imeas = a*d*(b*Isc + c)
    % Solve for Isc:
    %   Isc = Imeas/(a*d*b) - c/b = I/scale - offset
    %   where scale = a*d*b, offset = c/b
    %
    % I want to additionally multiply Isc by scaleFactor, which means I need to
    % do b -> b/scaleFactor.
    %
    % Also, I want to subtract Iincoh from Isc, which means I need to do...
    %
    %    Isc = Imeas/(a*d*b/scaleFactor) - c/(b/scaleFactor) - Incoh
    %
    %    offset_new = [c + Incoh*(b/scaleFactor)]/(b/scaleFactor)
    %    c -> c + Incoh*(b/scaleFactor);
    
    NewScalingModel = ScalingModel;
    
    for j=1:length(ScalingModel)
        SM = ScalingModel(j);
        
        [sgrid,ngrid] = ndgrid(SM.Ic.wx,1:length(SM.Ic.wy));
        
        Iincoh = 0;
        for n=1:numel(Molecules)
            [~,b] = Molecules(n).scatteringFactors(sgrid,'interp');
            Iincoh = Iincoh + b*occupancies(n);
        end
        
        NewScalingModel(j).c = SM.c + Iincoh.*SM.b(ngrid)/scaleFactor;
        NewScalingModel(j).b = SM.b/scaleFactor;
    end
    
    % END MAIN SCRIPT
    
    BatchProcess.saveToMatFile(options.scaleOut,...
        'options',options,'ScalingModel',NewScalingModel);
    
    BatchProcess.finish; % done!
    
catch errorMessage
    BatchProcess.stop(errorMessage);
end

tf = BatchProcess.hasCompleted;
errorMessage = BatchProcess.errorMessage;

end