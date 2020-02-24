function Wk = onePhononStructureFactors(Basis,Ops,hklGrid,Fgrid,Fxgrid,Fygrid,Fzgrid)

B = symm.AffineTransformation(Basis.orthogonalizationMatrix,[0;0;0]);
K1 = hklGrid(1).invert.invert; %bz grid

sGrid = arrayfun(@(g) latt.LatticeGrid(g,Basis.invert),hklGrid);

Wk = repmat({zeros(numel(sGrid),6*numel(Ops))},K1.N);

for j=1:numel(sGrid)

[sx,sy,sz] = sGrid(j).grid();

qx = 2*pi*sx;
qy = 2*pi*sy;
qz = 2*pi*sz;

Wgn = cell(1,numel(Ops));
for n=1:numel(Ops)

Op = B*Ops(n)*inv(B);

qr = [qx(:),qy(:),qz(:)]*Op.r;

F1 = Fgrid{j,n};
F1x = Fxgrid{j,n};
F1y = Fygrid{j,n};
F1z = Fzgrid{j,n};

% W = [ F*qx, F*qy, F*qz, Fy*qz - Fz*qy, Fz*qx - Fx*qz, Fx*qy - Fy*qx]
Wgn{1,n} = [F1(:).*qr(:,1),... % F*qx
    F1(:).*qr(:,2),...
    F1(:).*qr(:,3),...
    F1y(:).*qr(:,3) - F1z(:).*qr(:,2),... % Fy*qz - Fz*qy
    F1z(:).*qr(:,1) - F1x(:).*qr(:,3),... % Fz*qx - Fx*qz
    F1x(:).*qr(:,2) - F1y(:).*qr(:,1)]; % Fx*qy - Fy*qx
end
W = cell2mat(Wgn);

% calculate indices [n1,n2,n3] in bz array for each point in the input structure
[h1,h2,h3] = hklGrid(j).grid();
[n1,n2,n3] = K1.frac2ind(h1,h2,h3);
ix = sub2ind(K1.N,n1,n2,n3);

% assign wk
for k=1:numel(Wk)
    Wk{k}(j,:) = W(ix(k),:);
end

end