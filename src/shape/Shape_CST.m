%% Shape_CST Class
% This class defines the behavior of a constant strain triangular (CST)
% isoparametric element. It provides methods for evaluating shape functions,
% computing Jacobian matrices, derivatives, strain-displacement matrices, 
% and integration points for finite element analysis.
% 
%% Methods:
% * *shapeFnc*: Evaluate the shape functions at a given point in the 
%               natural coordinate system.
% * *shapeFncMtrx*: Get the shape function matrix.
% * *NuMtrx*: Compute the shape function matrix for displacement 
%             interpolation.
% * *shapeFncDrv*: Compute the derivatives of the shape functions with 
%                  respect to natural coordinates.
% * *JacobianMtrx*: Compute the Jacobian matrix for the element.
% * *detJacobian*: Compute the determinant of the Jacobian matrix.
% * *dNdxMatrix*: Compute the derivatives of the shape functions with 
%                 respect to global coordinates.
% * *BMatrix*: Compute the strain-displacement matrix.
% * *getIntegrationPointsCST*: Get integration points and weights for CST 
%                              elements based on quadrature order.
% * *getIntegrationPoints*: Get integration points and weights, with 
%                           optional subdivision for fractured elements.
% * *coordNaturalToCartesian*: Transform a point from natural to global 
%                              Cartesian coordinates.
% * *coordCartesianToNatural*: Transform a point from global Cartesian to 
%                              natural coordinates.
% * *areaTriangle*: Compute the area of a triangle given its vertices.
% * *getSizeGramMtrx*: Get the size of the Gram matrix.
% * *integrandGramMtrx*: Compute the integrand for the Gram matrix.
% * *getSizeStressIntVct*: Get the size of the stress interpolation vector.
% * *polynomialStress*: Define the polynomial stress interpolation.
% * *integrandStressIntVct*: Compute the integrand for the stress 
%                            interpolation vector.
% * *getTriangleLinearQuadrature*: Get integration points and weights for 
%                                  linear quadrature in a triangle.
% 
%% Author
% - Danilo Cavalcanti
%
%% Version History
% - Version 1.00: Initial version (January 2023).
%
%% Class Definition
classdef Shape_CST < Shape
    %% Constructor method
    methods
        function this = Shape_CST()
            this = this@Shape('CST');
        end
    end

    %% Public methods: methods defined in the abstract superclass
    methods
        % -----------------------------------------------------------------
        % Evaluate the shape function at a given point X of a linear 
        % quadrilateral isoparametric element.
         function N = shapeFnc(~,Xn)
            % Natural coordinates of the given point
            r = Xn(1); s = Xn(2);

            % Shape functions
            N1 = 1.0 - r - s;
            N2 = r;
            N3 = s;
            N   = [ N1  N2  N3 ];
         end

         % -----------------------------------------------------------------
         % Get the shape function matrix
         function Nm = shapeFncMtrx(this,Xn)
             % Vector with the shape functions
             Nm = this.shapeFnc(Xn);
         end

         % -----------------------------------------------------------------
         % Get the shape function matrix
         function Nu = NuMtrx(~,N)
             % Vector with the shape functions
             Nu = [N(1)  0.0   N(2)  0.0   N(3)  0.0;
                   0.0   N(1)  0.0   N(2)  0.0   N(3)];
         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions wrt to the
         % natural coordinate s
         function dNdxn = shapeFncDrv(~,~)
            % Derivatives of the shape functions
            dN1_dr = -1.0;    dN1_ds = -1.0;
            dN2_dr =  1.0;    dN2_ds =  0.0;
            dN3_dr =  0.0;    dN3_ds =  1.0;
            dNdxn   = [ dN1_dr  dN2_dr  dN3_dr;
                        dN1_ds  dN2_ds  dN3_ds];
         end

         % -----------------------------------------------------------------
         % Compute the jacobian matrix
         function J = JacobianMtrx(this,X,Xn)
            % Compute the shape function derivatives wrt to the natural
            % coordinate system
            dNdxn = this.shapeFncDrv(Xn);
            
            % Jacobian matrix
            J = dNdxn*X;
         end

         % -----------------------------------------------------------------
         % Compute the determinant of the jacobian
         function detJ = detJacobian(this,X,Xn)
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);
            detJ = det(J);
         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions matrix
         function [dNdx,detJ] = dNdxMatrix(this,X,Xn)
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);

            % Determinant of the Jacobian matrix
            detJ = det(J);

            % Compute the derivatives of the shape functions wrt to the
            % natural coordinate system
            dNdxn = this.shapeFncDrv(Xn);

            % Compute the derivatives of the shape functions wrt to the
            % global cartesian coordinate system
            dNdx = J\dNdxn;
         end

         % -----------------------------------------------------------------
         % Compute the strain-displacement matrix
         function [B] = BMatrix(~,dNdx)
            B = zeros(4,3*2);
            for i = 1:3
                B(1,2*i-1) = dNdx(1,i); 
                B(2,2*i)   = dNdx(2,i);
                B(4,2*i-1) = dNdx(2,i); B(4,2*i) = dNdx(1,i);
            end
         end

         % -----------------------------------------------------------------
         % Get the integration points:
         % Output:
         %      X: Coordinates of the integration points in the natural
         %         system
         %      w: Weight associated to each integration point
         % %      nIntPoints: Total number of integration points
         function [X,w,nIntPoints] = getIntegrationPointsCST(~,order)
             if order == 2
                 X          = [1/6,  1/6;
                               2/3,  1/6;
                               1/6,  2/3]';  
                 w          = [1/6,1/6,1/6];
                 nIntPoints = 3;
             elseif order == 1
                 X          = [1/3, 1/3]';
                 w          = 0.5;
                 nIntPoints = 1;
             end
         end

         %------------------------------------------------------------------
        % Function to get the matrices with the coordinates and weights of
        % integration points in the natural coordinate system. 
        % The definition of the integration points depends on the choice of
        % the order of the quadrature rule and the option for subdividing
        % the domain or not.
        % Output:
        %
        %   X : matrix with the coordinates of the integration points in
        %       the natural coordinate system. Each column of this matrix
        %       gives the coordinates of a point.
        %   W : vector with the weights associated with each point
        %   n : number of integration points
        %
        function [X,W,n] = getIntegrationPoints(this,intOrder,elem)
            if nargin == 2
                subDivInt = false;
            else
                subDivInt = elem.subDivInt;
            end

            if subDivInt == false
                [X,W,n] = this.getIntegrationPointsCST(intOrder);
            else
                elemNodes = [elem.node]; FractSeg = [];
                nnodes = size(elem.node,1);
                k = 1;
                for i = 1:elem.nfrac
                    % Number of fracture segments
                    nfracSeg = length(elem.fracture{i});
                    % Get the fracture tip nodes
                    fractiNodes = zeros(2);
                    fractiNodes(1,:) = elem.fracture{i}(1).node(1,:);
                    fractiNodes(2,:) = elem.fracture{i}(nfracSeg).node(2,:);
                    % Add to the element nodes
                    elemNodes = [elemNodes;fractiNodes];
                    FractSeg  = [FractSeg;nnodes+k,nnodes+k+1];
                    k = k + 2;
                end
                
                % Perform a Delaunay triangulation to create subelements
                DT = delaunayTriangulation(elemNodes,FractSeg);

                % Number of subelements
                nSubElem = size(DT.ConnectivityList,1);

                % Get the integration points of a triangle
                shapeCST = Shape_CST();
                [x,w,nIntPoints] = shapeCST.getIntegrationPointsCST(intOrder);

                % Total number of integration points
                n = nIntPoints * nSubElem;

                % Initialize the matrix of coordinates and the vector of
                % weights
                W = zeros(1,n);
                X = zeros(2,n);

                % Fill the matrix X and the vector W
                count = 1;
                for subelem = 1:nSubElem
                    for ip = 1:length(w)
                        % Find the global cartesian coordinates of the integration point
                        [Xi] = shapeCST.coordNaturalToCartesian(...
                            DT.Points(DT.ConnectivityList(subelem,:),:),x(:,ip));

                        % Weight
                        detJel = this.detJacobian(elem.node,x(:,ip));
                        detJ   = shapeCST.detJacobian(DT.Points(DT.ConnectivityList(subelem,:),:),x(:,ip));
                        W(count) = w(ip)*detJ/detJel;
                
                        % Mapping the integration point to the quadrilateral
                        % element natural system
                        X(:,count) = this.coordCartesianToNatural(elem.node,Xi)';

                        count = count + 1;

                    end
                end 
            end
        end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of an
         %           CST element
         %   Xn   : vector with the xi and eta coordinates of a point in 
         %          the natural coordinate system
         %
         % Output:
         %   X : vector with the x and y coordinates of a point in the 
         %       global coordinate system
         function X = coordNaturalToCartesian(this,NODE,Xn)
            % Extract the nodal coordinates
            x = NODE(:,1);
            y = NODE(:,2);
            
            % Vector with the shape functions
            Nv = this.shapeFnc(Xn);
            
            % Initialize output
            X = [0.0, 0.0];
            
            % Interpolation the position
            X(1) = Nv(1)*x(1) +  Nv(2)*x(2) +  Nv(3)*x(3);
            X(2) = Nv(1)*y(1) +  Nv(2)*y(2) +  Nv(3)*y(3);
         end

         % ----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Reference: Section 7.3 from the book Concepts and Applications
         % of Finite Element Analysis from Cook et al. (2001)
         function Xn = coordCartesianToNatural(this,NODE,X)
            % Compute the area of the triangle
            A = this.areaTriangle(NODE(1,:), NODE(2,:), NODE(3,:));

            % Compute the areas of the sub-triangles
            A1 = this.areaTriangle(NODE(2,:), NODE(3,:), X);
            A2 = this.areaTriangle(NODE(3,:), NODE(1,:), X);

            % Compute the area coordinates
            xi1 = A1/A;
            xi2 = A2/A;
            xi3 = 1.0 - xi1 - xi2;

            Xn = [xi2, xi3];
         end

        % -----------------------------------------------------------------
        % Compute the area of a triangle given the three point in
        % counterclock-wise order
        function A = areaTriangle(~, P1, P2, P3)
            P1P2 = P2 - P1;
            P1P3 = P3 - P1;
            A = (P1P2(1)*P1P3(2) - P1P2(2)*P1P3(1))*0.5;
        end

        % -----------------------------------------------------------------
        % Size of the Gram matrix
        % The stress field in a CST element is constant
        function n = getSizeGramMtrx(~)
            n = 1;
        end

        % -----------------------------------------------------------------
        % Integrand to compute the Gram Matrix
        % The definition of this matrix is associated to the order of the
        % stress field inside the element domain. For a linear
        % triangular isoparametric element, the stress field will be
        % assumed to have a constant
        function dH = integrandGramMtrx(~,~,~)
            dH = 1.0;    
        end

        % -----------------------------------------------------------------
        % Size of the stress interpolation vector
        function n = getSizeStressIntVct(~)
            n = 1;
        end

        % -----------------------------------------------------------------
        % Polynomial stress interpolation
        function p = polynomialStress(~,~)
            p = 1.0;
        end

        % -----------------------------------------------------------------
        % Integrand to compute the stress interpolation vector
        function dS = integrandStressIntVct(~,s,~,jumpOrder)
            if jumpOrder == 0
                dS = 1.0;
            elseif jumpOrder == 1
                dS = [ 1.0    s ];
            end
        end
    end

    methods(Static)
        function [X,w,nIntPoints] = getTriangleLinearQuadrature()
            % Coordinates of the integration points
            X          = [1/6,  1/6;
                          2/3,  1/6;
                          1/6,  2/3]';  
             % Weights
             w          = [1/6,1/6,1/6];
             % Total number of integration points
             nIntPoints = 3;
        end
    end
end
