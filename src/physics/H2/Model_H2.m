%% Model_H2 class
%
% Two-phase flow finite element model.
%
% Each node has two degrees of freedom (dof). The liquid phase pressure (p)
% and the gas phase pressure (pg).
%
%% Author
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class definition
classdef Model_H2 < Model_H
    %% Constructor method
    methods
        function this = Model_H2()
            this = this@Model_H(false);
            this.ndof_nd = 2;       % Number of dofs per node
            this.physics = 'H2';    % Tag with the physics name
            disp("*** Physics: Two-phase hydraulic (H2)");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function setMaterial(this,porousMedia,liquidFluid,gasFluid)
            if nargin < 4
                disp('Error in setMaterial: insuficient number of inputs.');
                disp('Physics H2 requires 3 attribute(s): porousMedia, liquidFluid, gasFluid.');
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
            this.mat = struct('porousMedia',porousMedia,'liquidFluid',liquidFluid,'gasFluid',gasFluid);
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
                pl_dofs = this.getElementDofs(el,1);
                pg_dofs = this.getElementDofs(el,2);
                elements(el) = RegularElement_H2(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder, pl_dofs, pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        function setGasPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 2, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 2, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 2, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 2, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 2, value);
        end

        % -----------------------------------------------------------------
        function setGasPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 2, value);
        end

        % -----------------------------------------------------------------
        function setInitialGasPressureAtDomain(this, value)
            this.setInitialDofAtDomain(2, value);
        end

        % -----------------------------------------------------------------
        function setInitialGasPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 2, value);
        end

        % -----------------------------------------------------------------
        function plotGasPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            FEMPlot = FEMPlot(this);
            FEMPlot.plotGasPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        function plotCapillaryPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            FEMPlot = FEMPlot(this);
            FEMPlot.plotCapillaryPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

    end

    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           Pl        Pg\n');
        end

    end
end