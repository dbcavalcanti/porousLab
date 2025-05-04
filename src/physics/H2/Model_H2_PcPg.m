%% Model_H2_PcPg Class
% This class extends the _Model_H2_ and represents a two-phase flow finite
% element model.
% Each node has two degrees of freedom:
% 
% * 1 capillary phase pressure (Pc)
% * 1 gas phase pressure (Pg)
%% Methods
% * *initializeElements*: Initializes the elements of the model with their 
%                         properties.
% * *setCapillaryPressureDirichletBCAtNode*: Sets pressure Dirichlet 
%                                            boundary conditions for 
%                                            gas-phase pore pressure at a 
%                                            specific node.
% * *setCapillaryPressureDirichletBCAtPoint*: Sets pressure Dirichlet 
%                                             boundary conditions for 
%                                             gas-phase pore pressure at a 
%                                             specific point.
% * *setCapillaryPressureDirichletBCAtBorder*: Sets pressure Dirichlet 
%                                              boundary conditions for 
%                                              gas-phase pore pressure at 
%                                              a specific border.
% * *setInitialCapillaryPressureAtDomain*: Sets the initial capillary 
%                                          pressure value across the 
%                                          entire domain.
% * *printResultsHeader*: Prints a header for the results table, showing
%                         node ID, and pressures (Pc, Pg).
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
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
        % Initializes the elements of the model with the corresponding
        % properties
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
                            this.NODE(this.ELEM{el},:), this.ELEM{el},...
                            this.t, emat, this.intOrder,pc_dofs,pg_dofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                if this.gravityOn
                    elements(el).type.gravityOn = true;
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        % Prescribe the capillary pressure Dirichlet boundary condition at 
        % a node
        function setCapillaryPressureDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, 1, value);
        end

        % -----------------------------------------------------------------
        % Prescribe the capillary pressure Dirichlet boundary condition at 
        % a point
        function setCapillaryPressureDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, 1, value);
        end

        % -----------------------------------------------------------------
        % Prescribe the capillary pressure Dirichlet boundary condition at 
        % a border
        function setCapillaryPressureDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, 1, value);
        end

        % -----------------------------------------------------------------
        % Sets the initial capillary pressure value for the whole domain
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