function Sol = solventMaskBlobs(Atoms,isIon,Basis,SpaceGroup,varargin)

    [xsol,ysol,zsol,LG] = latt.solventMaskPoints(Basis,SpaceGroup.generalPositions,...
        Atoms.x,Atoms.y,Atoms.z,Atoms.r,isIon,varargin{:});

    xyz = table2array(Atoms(:,{'x','y','z'}));
    N = size(Atoms,1);

    idxmin = knnsearch(xyz,[xsol,ysol,zsol]);
    ixlist = accumarray(idxmin,1:numel(idxmin),[N,1],@(d) {d});

    r1 = accumarray(idxmin,xsol,[N,1],@mean);
    r2 = accumarray(idxmin,ysol,[N,1],@mean);
    r3 = accumarray(idxmin,zsol,[N,1],@mean);

    v = det(LG.delta)*cellfun(@numel,ixlist);

    G = latt.GaussianDensitySum.empty();
    for n=1:numel(ixlist)

    un = [xsol(ixlist{n})-r1(n),ysol(ixlist{n})-r2(n),zsol(ixlist{n})-r3(n)];
    Un = un'*un/numel(ixlist{n});

    ellipsoidVolume = sqrt(det(Un)/.2.^3)*(4*pi/3);
    equivalentDensity = v(n)/ellipsoidVolume;

    G(n,1) = latt.GaussianDensitySum.approximateEllipsoid(Un,equivalentDensity); % occupancy = 1
    end

    Sol = table(r1,r2,r3,G,'VariableNames',{'xsol','ysol','zsol','fsol'});
end