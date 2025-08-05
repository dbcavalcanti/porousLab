%% regularMesh Function
% This function generates a regular finite element mesh with quadrilateral
% elements. It allows for optional quadratic distribution of nodes in the
% x and y directions and supports both linear quadrilateral (ISOQ4) and
% constant strain triangle (CST) element types.
% 
%% Inputs
% * *Lx*: Length in the x-direction (scalar)
% * *Ly*: Length in the y-direction (scalar)
% * *Nx*: Number of elements in the x-direction (integer)
% * *Ny*: Number of elements in the y-direction (integer)
% * *xo*: Fixed x-coordinates (optional, default: [])
% * *yo*: Fixed y-coordinates (optional, default: [])
% * *type*: Type of finite element ('ISOQ4' for quadrilateral or 'CST' for 
%           constant strain triangle, default: 'ISOQ4')
%
%% Outputs
% * *Node*: Matrix of node coordinates (Nx*Ny x 2), where each row 
%           contains the x and y coordinates of a node.
% * *ELEM*: Cell of element connectivity. Each row defines the nodes
%           forming an element. For 'ISOQ4', each row has 4 nodes; 
%           for 'CST', each row has 3 nodes.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Function definition
function [Node,ELEM] = regularMesh(Lx,Ly,Nx,Ny,xo,yo,type,cx,Ax,cy,Ay)

if nargin < 5, xo = []; yo = []; end
if nargin < 7, type = 'ISOQ4'; end
if nargin < 9, cx = 0; Ax = 0; end
if nargin < 11, cy = 0; Ay = 0; end

% Uniform parameter in [0,1]
s = linspace(0,1,Nx+1);  

% width of the dip
sigma = max(min(8/Nx,1.0),0.0);

% 2) inverted‐Gaussian “dip” PDF
pdf = 1 - Ax*exp( -((s - cx).^2) ./ (2*sigma^2) );

% 3) normalize PDF so area=1 (optional but cleaner)
pdf = pdf / trapz(s,pdf);

% 4) CIDF warp
ux = cumtrapz(s, pdf);
ux = ux / ux(end);     % ensure u(1)=1 exactly
xcoord = getUniquePoints(Lx*ux, xo, 0.2*Lx/Nx);

% Uniform parameter in [0,1]
s = linspace(0,1,Ny+1);   

% width of the dip
sigma = max(min(10/Ny,1.0),0.0);

% 2) inverted‐Gaussian “dip” PDF
pdf = 1 - Ay*exp( -((s - cy).^2) ./ (2*sigma^2) );

% 3) normalize PDF so area=1 (optional but cleaner)
pdf = pdf / trapz(s,pdf);

% 4) CIDF warp
uy = cumtrapz(s, pdf);
uy = uy / uy(end);     % ensure u(1)=1 exactly
ycoord = getUniquePoints(Ly*uy, yo, 0.2*Ly/Ny);

% Number of nodes in each direction
Nx = length(xcoord) - 1;
Ny = length(ycoord) - 1;

% Coordinates of the nodes
if Nx < Ny
    % Number the nodes first along x and then along y
    [Y,X]= meshgrid(ycoord,xcoord);
else
    % Number the nodes first along y and then along x
    [X,Y]= meshgrid(xcoord,ycoord);
end
Node= [reshape(X,numel(X),1) reshape(Y,numel(Y),1)];

% Initialize the element matrix
if strcmp(type,'CST')
    ELEM= cell(2*Nx*Ny,1);
else
    ELEM= cell(Nx*Ny,1);
end

% Numbering elements first along x and then along y
k= 1;
if Nx < Ny 
    for j=1:Ny
        for i=1:Nx  
            n1 = (j-1)*(Nx+1)+i; n2 = j*(Nx+1)+i;
            ELEM{k} = [n1, n1+1, n2+1, n2];
            k = k+1;
        end
    end
else
    for j=1:Ny
        for i=1:Nx  
            n1 = (i-1)*(Ny+1)+j; n2 = i*(Ny+1)+j;
            if strcmp(type,'CST')
                ELEM{k}   = [n1 n2 n2+1];
                ELEM{k+1} = [n2+1 n1+1 n1];
                k = k+2;
            else
                ELEM{k} = [n1 n2 n2+1 n1+1];
                k = k+1;
            end
        end
    end
end

end

% -------------------------------------------------------------------------
% Filters and combines unique points from two sets of points
function x = getUniquePoints(x0,xfix,tol)
if isempty(xfix) == true
    x = x0;
    return
end
x = xfix;
for i = 1:length(x0)
    dx = abs(x0(i) - xfix);
    if any(dx < tol) == false
        x = [x, x0(i)];
    end
end
x = sort(x);
end