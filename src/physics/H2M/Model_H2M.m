%% Model_H2M class
%
% Hydromechanical wiht two-phase fluid flow finite element model.
%
% Each node has 3 degrees of freedom (dof):
% - 2 displacement components (ux,uy)
% - 1 liquid-phase pore pressure (pl)
% - 1 gas-phase pore pressure (pg)
%
%% Author
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class definition
classdef Model_H2M < Model_M    
    %% Constructor method
    methods
        function this = Model_H2M()
            this = this@Model_M(false);
            this.ndof_nd = 4;        % Number of dofs per node
            this.physics = 'H2M';    % Tag with the physics name
            disp("*** Physics: Hydromechanical with two-phase flow");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function setMaterial(this,porousMedia,liquidFluid,gasFluid)
            if nargin < 4
                disp('Error in setMaterial: insuficient number of inputs.');
                disp('Physics H2M requires 3 attribute(s): porousMedia, liquidFluid, gasFluid.');
                error('Error in setMaterial.');
            end
            if ~isa(porousMedia,'PorousMedia')
                disp('Error in setMaterial: porousMedia is not a PorousMedia object.');
                error('Error in setMaterial.');
            end
            if ~isa(liquidFluid,'Fluid')
                disp('Error in setMaterial: liquidFluid is not a Fluid object.');
                error('Error in setMaterial.');
            end
            if ~isa(gasFluid,'Fluid')
                disp('Error in setMaterial: gasFluid is not a Fluid object.');
                error('Error in setMaterial.');
            end
            this.mat  = struct('porousMedia',porousMedia,'liquidFluid',liquidFluid,'gasFluid',gasFluid);
        end

        %------------------------------------------------------------------
        function initializeElements(this)
            % Initialize the vector with the Element's objects
            elements(this.nelem,1) = Element(); 

            % Assemble the properties to the elements' objects
            for el = 1 : this.nelem
                % Create the material for the element
                emat =struct( ...
                        'porousMedia',this.mat.porousMedia(this.matID(el)), ...
                        'liquidFluid',this.mat.liquidFluid,...
                        'gasFluid',this.mat.gasFluid);
                udofs = this.getElementDofs(el,[1,2]);
                pl_dofs = this.getElementDofs(el,3);
                pg_dofs = this.getElementDofs(el,4);
                elements(el) = RegularElement_H2M(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,udofs,pl_dofs,pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                            this.isPlaneStress);
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        function setPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 3, value);
        end

        % -----------------------------------------------------------------
        function setPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 3, value);
        end

        % -----------------------------------------------------------------
        function setPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 3, value);
        end

        % -----------------------------------------------------------------
        function setPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 3, value);
        end

        % -----------------------------------------------------------------
        function setPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 3, value);
        end

        % -----------------------------------------------------------------
        function setPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 3, value);
        end

        % -----------------------------------------------------------------
        function setInitialPressureAtDomain(this, value)
            this.setInitialDofAtDomain(3, value);
        end

        % -----------------------------------------------------------------
        function setInitialPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 3, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 4, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 4, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 4, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 4, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 4, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 4, value);
        end

        % -----------------------------------------------------------------
        function setInitialGasPressureAtDomain(this, value)
            this.setInitialDofAtDomain(4, value);
        end

        % -----------------------------------------------------------------
        function setInitialGasPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 4, value);
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotGasPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotGasPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotCapillaryPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotCapillaryPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

    end
        %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           ux        uy        Pl        Pg\n');
        end

    end
end