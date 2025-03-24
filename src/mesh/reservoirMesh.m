function [NODE] = reservoirMesh(NODE, XD, SEGD)

% Fault tangential vector
m = (XD(SEGD(2),:) - XD(SEGD(1),:))'/ norm(XD(SEGD(2),:) - XD(SEGD(1),:));

% Vertical vector
ey = [0.0; 1.0];

% Angle between the tangential vector and the vertical vector
theta = acos(dot(m,ey));

% Bottom coordinates of the fault
xf0 = XD(SEGD(1),1);

% x-coordinates of the regular mesh
xp = unique(NODE(:,1));

% Number of nodes in each direction
nx = length(xp);
ny = length(unique(NODE(:,2)));

% Number of elements in the x-direction
nelemx = nx - 1;

% Horizontal length of the domain 
Lx = max(NODE(:,1));

% Horizontal length of the elements in the regular mesh
lex = Lx/nelemx;

% x-coordinate of the center point of the fault
xfc = (XD(SEGD(2),1) + XD(SEGD(1),1))*0.5;

% Coordinates of the nodes
nleft  = length(find(xp < xfc));
nright = nx - nleft;

% Changing the x-coordinates of the nodes
for j = 1:ny

    % y-coordinate of the nodes
    yj = NODE(j,2);

    % x-coordinate of the fault at this height
    xf = xf0 + tan(theta)*yj;

    % Id of the nodes with this y-coordinate
    id = linspace(j,j+ny*(nx-1),nx);

    % New coordinates of the nodes on the left and on the right side of the
    % fault
    xleft  = linspace(0.0, xf-lex/2, nleft);
    xright = linspace(xf+lex/2, Lx, nright);
    xnew   = [xleft , xright];

    % Update the x-coordinates at the NODE matrix
    NODE(id,1) = xnew;
end

end