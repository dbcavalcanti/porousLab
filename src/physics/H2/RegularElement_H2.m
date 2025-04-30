%% RegularElement_H2 Class
% This class defines a finite element for a two-phase flow formulation 
% using the liquid pressure (Pl) and the gas pressure (Pg) as primary 
% variables. It extends the _RegularElement_ class and provides methods 
% for initializing integration points, assembling element matrices and 
% vectors, and computing various fields such as pressure and saturation.
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          element using the shape function and material 
%                          properties.
% * *elementData*: Assembles the element stiffness matrix, damping matrix, 
%                  internal force vector, external force vector, and 
%                  derivative of internal force with respect to 
%                  displacement.
% * *permeabilityTensors*: Computes the permeability tensors for the 
%                          element.
% * *compressibilityCoeffs*: Computes the compressibility coefficients 
%                            for the element.
% * *lumpedCompressibilityMatrix*: Computes the lumped compressibility 
%                                  matrices based on the element volume 
%                                  and compressibility coefficients.
% * *addGravityForces*: Adds the contribution of gravity forces to the 
%                       external force vector.
% * *getNodalLiquidPressure*: Retrieves the nodal liquid pressure values.
% * *getNodalGasPressure*: Retrieves the nodal gas pressure values.
% * *getNodalCapillaryPressure*: Retrieves the nodal capillary pressure 
%                                values.
% * *pressureField*: Computes the pressure field at a given position 
%                    inside the element.
% * *gasPressureField*: Computes the gas pressure field at a given 
%                       position inside the element.
% * *capillaryPressureField*: Computes the capillary pressure field at a 
%                             given position inside the element.
% * *liquidSaturationField*: Computes the liquid saturation field at a 
%                            given position inside the element.
% * *gasSaturationField*: Computes the gas saturation field at a given 
%                         position inside the element.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef RegularElement_H2 < RegularElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glp        = [];            
        glpg       = [];            % Vector of the regular degrees of freedom
        nglp       = 0;             % Number of regular p-dof
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_H2(node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement(node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.glp      = glp;
            this.glpg     = glpg;
            this.gle      = [glp , glpg];
            if (length(this.glp) ~= length(this.glpg))
                error('Wrong number of pressure dofs');
            end
            this.nglp     = length(this.glp);
            this.ngle     = length(this.gle);
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(this.intOrder);

            % Initialize the integration points objects
            
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = Material_H2(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        %
        % Output:
        %    Ke : element "stiffness" matrix
        %    Ce : element "damping" matrix
        %    fe : element "external force" vector
        %    fi : element "internal force" vector
        % dfidu : element matrix of derivative of the internal force with 
        %         respect to displacement
        %
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

            % Initialize the sub-matrices
            Hll = zeros(this.nglp, this.nglp);
            Hlg = zeros(this.nglp, this.nglp);
            Hgl = zeros(this.nglp, this.nglp);
            Hgg = zeros(this.nglp, this.nglp);  
            Sll = zeros(this.nglp, this.nglp);
            Slg = zeros(this.nglp, this.nglp);
            Sgl = zeros(this.nglp, this.nglp);
            Sgg = zeros(this.nglp, this.nglp);
            dfidu = zeros(2*this.nglp, 2*this.nglp);

            % Initialize external force vector
            fel = zeros(this.nglp, 1);
            feg = zeros(this.nglp, 1);

            % Initialize internal force vector
            fi = zeros(2*this.nglp, 1);
            
            % Vector of the nodal pore-pressure dofs
            pc = this.getNodalCapillaryPressure();
            pg = this.getNodalGasPressure();

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Pressure values at the integration point
                pcIP = Np * pc;
                pgIP = Np * pg;

                % Compute the saturation degree at the integration point
                Sl = this.intPoint(i).constitutiveMdl.saturationDegree(pcIP);
        
                % Compute the permeability matrix
                [kll, klg, kgl, kgg] = this.permeabilityTensors(this.intPoint(i),pgIP,pcIP,Sl);

                % Get compressibility coefficients
                [cll, clg, cgl, cgg] = this.compressibilityCoeffs(this.intPoint(i),pgIP,pcIP,Sl);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
        
                % Compute permeability sub-matrices
                Hll = Hll + Bp' * kll * Bp * c;
                Hlg = Hlg + Bp' * klg * Bp * c;
                Hgl = Hgl + Bp' * kgl * Bp * c;
                Hgg = Hgg + Bp' * kgg * Bp * c;

                % Compute compressibility matrices
                if ((this.massLumping) && (this.lumpStrategy == 1))
                    Sll = Sll + diag(cll*Np*c);
                    Sgg = Sgg + diag(cgg*Np*c);
                    Slg = Slg + diag(clg*Np*c);
                    Sgl = Sgl + diag(cgl*Np*c);
                elseif (this.massLumping == false)
                    Sll = Sll + Np' * cll * Np * c;
                    Sgg = Sgg + Np' * cgg * Np * c;
                    Slg = Slg + Np' * clg * Np * c;
                    Sgl = Sgl + Np' * cgl * Np * c;
                end
                
                % Compute the gravity forces
                if (this.gravityOn)
                    [fel,feg] = this.addGravityForces(fel,feg,Bp,kll,kgg,pgIP-pcIP,pgIP,c);
                end

                % Compute the element volume
                vol = vol + c;
            end

            % Compute the lumped mass matrix
            if ((this.massLumping) && (this.lumpStrategy == 2))
                [Sll,Slg,Sgl,Sgg] = lumpedCompressibilityMatrices(this, pc, pg, vol);
            end

            % Assemble the element permeability
            Ke = [ Hll, Hlg;
                   Hgl, Hgg ];

            % Assemble the element compressibility matrix
            Ce = [ Sll , Slg;
                   Sgl , Sgg ];

            % Assemble element external force vector
            fe = [fel; feg];
            
        end

        % -----------------------------------------------------------------
        % Compute the permeability tensors
        function [kll, klg, kgl, kgg] = permeabilityTensors(~,ip,pg,pc,Sl)
             [kll, klg, kgl, kgg] = ip.constitutiveMdl.permeabilityMtrcs(Sl,pg-pc,pg);
        end

        % -----------------------------------------------------------------
        % Compute the compressibility coefficients
        function [cll, clg, cgl, cgg] = compressibilityCoeffs(~,ip,pg,pc,Sl)
             [cll, clg, cgl, cgg] =  ip.constitutiveMdl.compressibilityCoeffs(Sl,pg-pc,pg);
        end

        %------------------------------------------------------------------
        % Compute the lumped mass matrices
        function [Sll,Slg,Sgl,Sgg] = lumpedCompressibilityMatrices(this, pc, pg, vol)

            % Shape function matrix
            Np = this.shape.shapeFncMtrx([0.0,0.0]);

            % Pressure values at the integration point
            pcIP = Np * pc;
            pgIP = Np * pg;

            % Compute the saturation degree at the integration point
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pcIP);

            % Get compressibility coefficients
            [cll, clg, cgl, cgg] = this.compressibilityCoeffs(this.intPoint(1),pgIP,pcIP,Sl);

            % Mass distribution factor
            factor = vol / this.nnd_el;

            % Compressibility matrices
            Sll = cll * factor * eye(this.nglp,this.nglp);
            Slg = clg * factor * eye(this.nglp,this.nglp);
            Sgl = cgl * factor * eye(this.nglp,this.nglp);
            Sgg = cgg * factor * eye(this.nglp,this.nglp);

        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vct
        function [fel,feg] = addGravityForces(this,fel,feg,Bp,kl,kg,pl,pg,c)

            % Get gravity vector
            grav = this.g * this.mat.porousMedia.b;

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity();
            rhog = this.mat.gasFluid.getDensity();

            % Compute the contribution of the gravitational forces
            fel = fel + Bp' * kl * rhol * grav * c;
            feg = feg + Bp' * kg * rhog * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalLiquidPressure(this)
            pl = this.ue(1:this.nglp);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the gas pressure
        function pg = getNodalGasPressure(this)
            pg = this.ue(1+this.nglp:end);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the capillary pressure
        function pc = getNodalCapillaryPressure(this)
            pl = this.getNodalLiquidPressure();
            pg = this.getNodalGasPressure();
            pc = pg - pl;
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = pressureField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   p   : pressure evaluated in "X"
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Get nodal pressures
            pl = this.getNodalLiquidPressure();

            % capillary field
            p = Nm*pl;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = gasPressureField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   p   : pressure evaluated in "X"
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Get nodal pressures
            pg = this.getNodalGasPressure();

            % capillary field
            p = Nm*pg;
        
        end

        %------------------------------------------------------------------
        % Function to compute the pressure field inside a given element
        function p = capillaryPressureField(this,X,ue)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   p   : pressure evaluated in "X"
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Get nodal pressures
            pc = this.getNodalCapillaryPressure();

            % capillary field
            p = Nm*pc;
        
        end

        %------------------------------------------------------------------
        % Function to compute the liquid saturation field inside a given element
        function Sl = liquidSaturationField(this,X,ue)

            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Capillary pressure at the given point
            pc = Nm*this.getNodalCapillaryPressure();

            % Compute the liquid saturation degree
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pc);
        
        end

        %------------------------------------------------------------------
        % Function to compute the gas saturation field inside a given element
        function Sg = gasSaturationField(this,X,ue)
            
            if nargin > 2, this.ue = ue; end
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Capillary pressure at the given point
            pc = Nm*this.getNodalCapillaryPressure();

            % Compute the liquid saturation degree
            Sl = this.intPoint(1).constitutiveMdl.saturationDegree(pc);

            % Gas saturation degree
            Sg = 1.0 - Sl;
        
        end
    end
end