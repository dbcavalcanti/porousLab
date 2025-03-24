function [NODE] = reservoirMeshY(NODE, XD, SEGD)

% Fault tangential vector
m = (XD(SEGD(2),:) - XD(SEGD(1),:))'/ norm(XD(SEGD(2),:) - XD(SEGD(1),:));

% Vertical vector
ex = [1.0; 0.0];

% Angle between the tangential vector and the horizontal vector
theta = acos(dot(m,ex));

% Bottom coordinates of the fault
yf0 = XD(SEGD(1),2);

% x and y-coordinates of the regular mesh
xp = unique(NODE(:,1));
yp = unique(NODE(:,2));

% Number of nodes in each direction
nx = length(xp);
ny = length(yp);

% Number of elements in the y-direction
nelemy = ny - 1;

% Horizontal length of the domain 
Ly = max(NODE(:,2));

% Horizontal length of the elements in the regular mesh
ley = Ly/nelemy;

% y-coordinate of the center point of the fault
yfc = (XD(SEGD(2),2) + XD(SEGD(1),2))*0.5;

% Coordinates of the nodes
nbottom  = length(find(yp < yfc));
ntop = ny - nbottom;

% Changing the x-coordinates of the nodes
for j = 1:nx

    % y-coordinate of the nodes
    xj = xp(j);

    % x-coordinate of the fault at this height
    yf = yf0 + tan(theta)*xj;

    % Id of the nodes with this x-coordinate
    id = linspace((j-1)*ny+1,j*ny,ny);
    % id = (j-1)*ny+1:j*ny;

    % New coordinates of the nodes on the left and on the right side of the
    % fault
    ybottom  = linspace(0.0, yf-ley/2, nbottom);
    ytop = linspace(yf+ley/2, Ly, ntop);
    ynew   = [ybottom , ytop];

    % Update the x-coordinates at the NODE matrix
    NODE(id,2) = ynew;
end

end