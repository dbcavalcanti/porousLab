function [Node,ELEM] = regularMeshY(Lx,Ly,Nx,Ny,xo,yo,type,quadDistrX,quadDistrY)
% -------------------------------------------------------------------------
% This function generates a regular finite element mesh with quadrilateral
% elements.
%
% Input:
%   Lx:   length in the x-direction
%   Ly:   lenght in the y-direction
%   Nx:   number of elements in the x-direction
%   Ny:   number of elements in the y-direction
%   type: type of finite element (default: linear quadrilateral)
%
% Output:
%   NODE: matrix with the nodes coordinates
%   ELEM: matrix with the elements connectivity
% -------------------------------------------------------------------------

if nargin < 5, xo = []; yo = []; end
if nargin < 7, type = 'Q4'; end
if nargin < 9, quadDistrX = false; quadDistrY = false; end

% Generate the NODE matrix
if quadDistrX == false
    xcoord = getUniquePoints(linspace(0,Lx,Nx+1),xo, 0.2*Lx/Nx);
else
    xcoord = getUniquePoints(linspace(0,1,Nx+1),xo, 0.2*Lx/Nx);
    xcoord = xcoord.^2;
    xcoord = xcoord * Lx;
end
if quadDistrY == false
    ycoord = getUniquePoints(linspace(0,Ly,Ny+1),yo, 0.2*Ly/Ny);
else
    ycoord = getUniquePoints(linspace(0,1,Ny+1),yo, 0.2*Ly/Ny);
    ycoord = ycoord.^2;
    ycoord = ycoord * Ly;
end

Nx = length(xcoord) - 1;
Ny = length(ycoord) - 1;
[Y,X]= meshgrid(ycoord,xcoord);
Node= [reshape(X,numel(X),1) reshape(Y,numel(Y),1)];

% Generate the elements 
k= 1;
if strcmp(type,'CST')
    ELEM= zeros(2*Nx*Ny,3);
else
    ELEM= zeros(Nx*Ny,4);
end
% Numbering elements first along x and then along y
for j=1:Ny
    for i=1:Nx  
        n1 = (j-1)*(Nx+1)+i; n2 = j*(Nx+1)+i;
        ELEM(k,:) = [n1, n1+1, n2+1, n2];
        k = k+1;
    end
end

end

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