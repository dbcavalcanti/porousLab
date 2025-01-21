%% RegularElement_H2 class
%
% This class defines a finite element of a two-phase flow formulation using
% the liquid pressure (Pl) and the gas pressure (Pg) as primary
% variables.
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef RegularElement_H2M < RegularElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        glu        = [];            % Displacement dofs
        glp        = [];            % Liquid phase pressure dofs
        glpg       = [];            % Gas phase pressure dofs
        nglu       = 0;             % Number of regular u-dof
        nglp       = 0;             % Number of regular p-dof
        anm        = 'PlaneStrain'; % Analysis model
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_H2M(type, node, elem, t, ...
                mat, intOrder, glu, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement(type, node, elem, t, ...
                mat, intOrder, massLumping, lumpStrategy, ...
                isAxisSymmetric);
            this.glu      = glu;
            this.glp      = glp;
            this.glpg     = glpg;
            this.gle      = [glu, glp, glpg];
            if (length(this.glp) ~= length(this.glpg))
                error('Wrong number of pressure dofs');
            end
            this.nglu     = length(this.glu);
            this.nglp     = length(this.glp);
            this.ngle     = length(this.gle);
            if isPlaneStress
                this.anm = 'PlaneStress';
            end
            if isAxisSymmetric
                this.anm = 'AxisSymmetrical';
            end
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
                constModel = Material_H2M(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
                intPts(i).initializeMechanicalAnalysisModel(this.anm);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        %
        % Output:
        %   Ke : element "stiffness" matrix
        %   Ce : element "damping" matrix
        %   fe : element "internal force" vector
        %
        function [Ke, Ce, fi, fe] = elementData(this)

            % Initialize the sub-matrices
            Kuu = zeros(this.nglu, this.nglu);
            Hll = zeros(this.nglp, this.nglp);
            Hlg = zeros(this.nglp, this.nglp);
            Hgl = zeros(this.nglp, this.nglp);
            Hgg = zeros(this.nglp, this.nglp);
            Qul = zeros(this.nglp, this.nglu);
            Qug = zeros(this.nglp, this.nglu);
            Sll = zeros(this.nglp, this.nglp);
            Slg = zeros(this.nglp, this.nglp);
            Sgl = zeros(this.nglp, this.nglp);
            Sgg = zeros(this.nglp, this.nglp);
            Opu = zeros(this.nglp, this.nglu);
            Ouu = zeros(this.nglu, this.nglu);

            % Initialize external force vector
            feu = zeros(this.nglu, 1);
            fel = zeros(this.nglp, 1);
            feg = zeros(this.nglp, 1);

            % Initialize the internal force vector
            fiu = zeros(this.nglu, 1);
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();
            pc = this.getNodalCapillaryPressure();
            pg = this.getNodalGasPressure();
            pl = this.getNodalLiquidPressure();

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 0.0];

            % Initialize the volume of the element
            vol = 0.0;

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Np = this.shape.shapeFncMtrx(this.intPoint(i).X);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(Bp);

                % Pressure values at the integration point
                pcIP = Np * pc;
                pgIP = Np * pg;

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();

                % Compute the saturation degree at the integration point
                Sl = this.intPoint(i).constitutiveMdl.saturationDegree(pcIP);
        
                % Compute the permeability matrix
                [kll, klg, kgl, kgg] = this.permeabilityTensors(this.intPoint(i),pgIP,pcIP,Sl);

                % Get compressibility coefficients
                [cul, cug, cll, clg, cgl, cgg] = this.compressibilityCoeffs(this.intPoint(i),pgIP,pcIP,Sl);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end
                
                % Compute the stiffness sub-matrix
                Kuu = Kuu + Bu' * Duu * Bu * c;

                % Internal force vector
                fiu = fiu + Bu' * stress * c;

                % Compute the hydromechanical coupling matrices
                Qul = Qul + Np' * cul * m' * Bu * c;
                Qug = Qug + Np' * cug * m' * Bu * c;
        
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
                if (this.mat.porousMedia.gravityOn)
                    [feu,fel,feg] = this.addGravityForces(feu,fel,feg,Np,Bp,kll,kgg,pgIP-pcIP,pgIP,c);
                end

                % Compute the element volume
                vol = vol + c;
            end

            % Compute the lumped mass matrix
            if ((this.massLumping) && (this.lumpStrategy == 2))
                [Sll,Slg,Sgl,Sgg] = lumpedCompressibilityMatrices(this, pc, pg, vol);
            end

            % Assemble the element permeability
            Ke = [ Kuu , -Qul', -Qug'; 
                   Opu ,  Hll , Hlg;
                   Opu ,  Hgl , Hgg ];

            % Assemble the element compressibility matrix
            Ce = [ Ouu , Opu', Opu';
                   Qul , Sll , Slg ;
                   Qug , Sgl , Sgg ];

            % Add contribution of the pressure to the internal force vector
            fiu = fiu - Qul' * pl - Qug' * pg;
            fil = Hll * pl + Hlg * pg;
            fig = Hgl * pl + Hgg * pg;

            % Assemble element internal force vector
            fi = [fiu; fil; fig];

            % Assemble element external force vector
            fe = [feu; fel; feg];
            
        end

        % -----------------------------------------------------------------
        % Compute the permeability tensors
        function [kll, klg, kgl, kgg] = permeabilityTensors(~,ip,pg,pc,Sl)
             [kll, klg, kgl, kgg] = ip.constitutiveMdl.permeabilityMtrcs(Sl,pg-pc,pg);
        end

        % -----------------------------------------------------------------
        % Compute the compressibility coefficients
        function [cul, cug, cll, clg, cgl, cgg] = compressibilityCoeffs(~,ip,pg,pc,Sl)
             [cll, clg, cgl, cgg] =  ip.constitutiveMdl.compressibilityCoeffs(Sl,pg-pc,pg);
             [cul, cug] = ip.constitutiveMdl.mechanicalCompressibilityCoeffs(Sl);
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
        function [feu,fel,feg] = addGravityForces(this,feu,fel,feg,Np,Bp,kl,kg,pl,pg,c)

            % Get gravity vector
            grav = this.mat.porousMedia.g * this.mat.porousMedia.b;

            % Shape function matrix
            Nu = this.shape.NuMtrx(Np);

            % Get the porous matrix density
            rhos = this.mat.porousMedia.getDensity();

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity(pl);
            rhog = this.mat.gasFluid.getDensity(pg);

            % Compute the contribution of the gravitational forces
            feu = feu + Nu' * rhos * grav * c;
            fel = fel + Bp' * kl * rhol * grav * c;
            feg = feg + Bp' * kg * rhog * grav * c;
            
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the displacement
        function u = getNodalDisplacement(this)
            u = this.ue(1:this.nglu);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalLiquidPressure(this)
            a = this.nglu + 1;
            b = this.nglu + this.nglp;
            pl = this.ue(a:b);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the gas pressure
        function pg = getNodalGasPressure(this)
            a = this.nglu + this.nglp + 1;
            pg = this.ue(a:end);
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