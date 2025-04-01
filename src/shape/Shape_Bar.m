%% Shape_Bar Class
% This is class defines the behavior of a linear bar isoparametric element.
% It provides methods for evaluation the shape functions, computing
% derivatived, Jacobian matrices, stain-displacement matrices, and
% performing coordinate transformations between natural and Cartesian
% coordinate systems. Additionally, it supports integration point
% generation and operations related to stress interpolation and Gram matrix
% computation.
%
%% Methods
% This class provides the following methods:
% 
% * *shapeFnc*: Evaluates the shape functions at a given point in the 
%               natural coordinate system.
% * *shapeFncMtrx*: Computes the shape function matrix.
% * *shapeFncDrv*: Computes the derivatives of the shape functions with 
%                  respect to the natural coordinate.
% * *JacobianMtrx*: Computes the Jacobian matrix for the element.
% * *detJacobian*: Computes the determinant of the Jacobian matrix.
% * *dNdxMatrix*: Computes the derivatives of the shape functions with 
%                 respect to the global coordinate system and the 
%                 determinant of the Jacobian.
% * *BMatrix*: Computes the strain-displacement matrix.
% * *getIntegrationPoints*: Provides the integration points, weights, 
%                           and the total number of integration points for 
%                           a given integration order.
% * *coordNaturalToCartesian*: Transforms a point from the natural 
%                              coordinate system to the global Cartesian 
%                              coordinate system.
% * *coordCartesianToNatural*: Placeholder method to transform a point from 
%                              the global Cartesian coordinate system to 
%                              the natural coordinate system (not 
%                              implemented).
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (January 2023).
%
%% Class Definition
classdef Shape_Bar < Shape
    %% Constructor method
    methods
        function this = Shape_Bar()
            this = this@Shape('Bar');
        end
    end

    %% Public methods: methods defined in the abstract superclass
    methods
        % -----------------------------------------------------------------
        % Evaluate the shape function at a given point X of a linear 
        % quadrilateral isoparametric element.
         function N = shapeFnc(~,Xn)
            % Natural coordinates of the given point
            s = Xn(1);

            % Shape functions
            N1 = (1.0 - s)/2.0;
            N2 = (1.0 + s)/2.0;
            N   = [ N1  N2 ];
         end

         % -----------------------------------------------------------------
         % Get the shape function matrix
         function Nm = shapeFncMtrx(this,Xn)
             % Vector with the shape functions
             N = this.shapeFnc(Xn);
            
             % Shape function matrix
             Nm = [N(1)  0.0   N(2)  0.0 ;
                   0.0   N(1)  0.0   N(2)];
         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions wrt to the
         % natural coordinate s
         function dNdxn = shapeFncDrv(~,~)
            dNdxn   = [ -0.5  0.5];
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
            detJ = norm(J);
         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions matrix
         function [dNdx,detJ] = dNdxMatrix(this,X,Xn)
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);

            % Determinant of the Jacobian matrix
            detJ = norm(J);

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
            B = zeros(4,4*2);
            for i = 1:4
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
         %      nIntPoints: Total number of integration points
         function [X,w,nIntPoints] = getIntegrationPoints(~,intOrder,~)
             if intOrder == 0
                 X = 0.0;
                 w = 2.0;
                 nIntPoints = 1;
             elseif intOrder == 1
                 % X          = [-0.577350269189626 0.577350269189626];
                 X          = [-1.0, 1.0];  
                 w          = [ 1.0, 1.0];
                 nIntPoints = 2;
             elseif intOrder == 2
                 X          = [-1.0 , 0.0, 1.0];  
                 w          = [ 0.5, 1.0, 0.5];
                 nIntPoints = 3;
             elseif intOrder == 3
                 X          = [-1.0 , -1/3, 1/3, 1.0];  
                 w          = [ 1/3, 4/3 , 4/3, 1/3];
                 nIntPoints = 4;
             end
         end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of a
         %          bar element
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
            X(1) = Nv(1)*x(1) +  Nv(2)*x(2);
            X(2) = Nv(1)*y(1) +  Nv(2)*y(2);
         end

         % ----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % MUST BE IMPLEMENTED.
         function Xn = coordCartesianToNatural(~,~,~)
            Xn = [];  
         end
    end
end
