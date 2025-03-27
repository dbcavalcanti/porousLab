%% Model_H2 class
%
% Two-phase flow finite element model
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model_H2_PcPg < Model_H2    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Model_H2_PcPg()
            this = this@Model_H2();
            this.physics = 'H2_PcPg';
        end
    end
    
    %% Public methods
    methods

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
                pc_dofs = this.getElementDofs(el,1);
                pg_dofs = this.getElementDofs(el,2);
                elements(el) = RegularElement_H2_PcPg(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,pc_dofs,pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        function setCapillaryPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        function setCapillaryPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 1, value);
        end

        % -----------------------------------------------------------------
        function setCapillaryPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 1, value);
        end

        % -----------------------------------------------------------------
        function setInitialCapillaryPressureAtDomain(this, value)
            this.setInitialDofAtDomain(1, value);
        end

    end

    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           Pc        Pg\n');
        end

    end
end