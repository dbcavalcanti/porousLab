function [SUPP, LOAD, PRESCDISPL] = boundaryConditionsDisplacement(NODE,CoordSupp,CoordLoad,CoordPresc,Lx, Ly, Nx, Ny)

SUPP       = zeros(size(NODE,1),2);
LOAD       = zeros(size(NODE,1),2);
PRESCDISPL = zeros(size(NODE,1),2);

% --- Generate loads
epsx= Lx/Nx/10; epsy= Ly/Ny/10;
for i=1:size(CoordLoad,1)
    Px = CoordLoad(i,1);
    Py = CoordLoad(i,2);
    cx = CoordLoad(i,3);
    cy = CoordLoad(i,4);
    if cx<0         % If cx<0, generate load along x at nodes with y==cy
        k= find(abs(NODE(:,2)-cy)<=epsy);
        k0 = find((abs(NODE(:,1)-0.0)<=epsx) .* (abs(NODE(:,2)-cy)<=epsy));
        kf = find((abs(NODE(:,1)-Lx)<=epsx) .* (abs(NODE(:,2)-cy)<=epsy));
    elseif cy<0      % If cy<0, generate load along y at nodes with x==cx
        k= find(abs(NODE(:,1)-cx)<=epsx);
        k0 = find((abs(NODE(:,2)-0.0)<=epsx) .* (abs(NODE(:,1)-cx)<=epsy));
        kf = find((abs(NODE(:,2)-Ly)<=epsx) .* (abs(NODE(:,1)-cx)<=epsy));
    else
        k= find(abs(NODE(:,1)-cx)<=epsx & abs(NODE(:,2)-cy)<=epsy);
        k0 = -1;
        kf = -1;
    end
    for j= 1:length(k)
        if (k(j) == k0) || (k(j) == kf)
            LOAD(k(j),:)= [Px Py]*0.5;
        else
            LOAD(k(j),:)= [Px Py];
        end
    end
end

% --- Generate prescribe pressure
epsx= Lx/Nx/10; epsy= Ly/Ny/10;
for i=1:size(CoordPresc,1)
    Px = CoordPresc(i,1);
    Py = CoordPresc(i,2);
    cx = CoordPresc(i,3);
    cy = CoordPresc(i,4);
    if cx<0         % If cx<0, generate load along x at nodes with y==cy
        k= find(abs(NODE(:,2)-cy)<=epsy);
    elseif cy<0      % If cy<0, generate load along y at nodes with x==cx
        k= find(abs(NODE(:,1)-cx)<=epsx);
    else
        k= find(abs(NODE(:,1)-cx)<=epsx & abs(NODE(:,2)-cy)<=epsy);
    end
    for j= 1:length(k)
        PRESCDISPL(k(j),:)= [Px Py];
    end
end

% --- Generate supports
for i=1:size(CoordSupp,1)
    rtx = CoordSupp(i,1); %if rtx==0,rtx= NaN; end
    rty = CoordSupp(i,2); %if rty==0,rty= NaN; end
    cx  = CoordSupp(i,3);
    cy  = CoordSupp(i,4);
    if cx<0          % If cx<0, generate load along x at nodes with y==cy
        k= find(abs(NODE(:,2)-cy)<=epsy);
    elseif cy<0      % If cy<0, generate load along y at nodes with x==cx
        k= find(abs(NODE(:,1)-cx)<=epsx);
    else
        k= find(abs(NODE(:,1)-cx)<=epsx & abs(NODE(:,2)-cy)<=epsy);
    end
    for j= 1:length(k)
        SUPP(k(j),:)= [rtx rty];
    end
end