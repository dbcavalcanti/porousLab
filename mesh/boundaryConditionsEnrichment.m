function LOAD_ENR_MID = boundaryConditionsEnrichment(NODE_D,CoordDischargeFract, Lx, Ly, Nx, Ny)

LOAD_ENR_MID = zeros(size(NODE_D,1),1);

% --- Generate loads
epsx= Lx/Nx/10; epsy= Ly/Ny/10;
for i=1:size(CoordDischargeFract,1)
    Px = CoordDischargeFract(i,1);
    cx = CoordDischargeFract(i,2);
    cy = CoordDischargeFract(i,3);
    if cx<0         % If cx<0, generate load along x at nodes with y==cy
        k= find(abs(NODE_D(:,2)-cy)<=epsy);
    elseif cy<0      % If cy<0, generate load along y at nodes with x==cx
        k= find(abs(NODE_D(:,1)-cx)<=epsx);
    else
        k= find(abs(NODE_D(:,1)-cx)<=epsx & abs(NODE_D(:,2)-cy)<=epsy);
    end
    for j= 1:length(k)
        LOAD_ENR_MID(k(j),:)= Px;
    end
end