%% Model_H class
%
% Single-phase fluid flow finite element model.
%
% Each node has one degree of freedom (dof) which is the pore-pressure (p)
%
%% Author
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class definition
classdef Model_H < Model    
    %% Constructor method
    methods
        function this = Model_H(printFlag)
            if nargin == 0, printFlag = true; end
            this = this@Model();
            this.ndof_nd = 1;       % Number of dofs per node
            this.physics = 'H';     % Tag with the physics name
            if (printFlag)
                disp("*** Physics: Single-phase hydraulic (H)");
            end
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function setMaterial(this,porousMedia,fluid)
            if nargin < 3
                disp('Error in setMaterial: insuficient number of inputs.');
                disp('Physics H requires 2 attribute(s): porousMedia, fluid.');
                error('Error in setMaterial.');
            end
            if ~isa(porousMedia,'PorousMedia')
                disp('Error in setMaterial: porousMedia is not a PorousMedia object.');
                error('Error in setMaterial.');
            end
            if ~isa(fluid,'Fluid')
                disp('Error in setMaterial: fluid is not a Fluid object.');
                error('Error in setMaterial.');
            end
            this.mat = struct('porousMedia',porousMedia,'fluid',fluid);
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
                        'fluid',this.mat.fluid);
                dof_e = this.getElementDofs(el,1);
                if (this.enriched == false)
                    elements(el) = RegularElement_H(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                else
                    elements(el) = EnrichedElement_H(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                end
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        function setPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        function setPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 1, value);
        end

        % -----------------------------------------------------------------
        function setPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 1, value);
        end

        % -----------------------------------------------------------------
        function setPressureNeumannBCAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        function setPressureNeumannBCAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, 1, value);
        end

        % -----------------------------------------------------------------
        function setPressureNeumannBCAtBorder(this, border, value)
            this.setNeumannBCAtBorder(border, 1, value);
        end

        % -----------------------------------------------------------------
        function setInitialPressureAtDomain(this, value)
            this.setInitialDofAtDomain(1, value);
        end

        % -----------------------------------------------------------------
        function setInitialPressureAtNode(this, nodeId, value)
            this.setInitialDofAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        function seg = initializeDiscontinuitySegArray(~,n)
            seg(n,1) = DiscontinuityElement_H([],[]);
        end

        % -----------------------------------------------------------------
        function seg = initializeDiscontinuitySegment(~,nodeD,matD)
            seg = DiscontinuityElement_H(nodeD,matD);
        end

        %------------------------------------------------------------------
        function mat = createMaterialDataStructure(this)
            mat = struct( ...
                'fluid',this.fluid, ...
                'initialAperture',this.initialAperture);
        end
    end

    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           Pl\n');
        end

    end
end