%% Model_HM class
%
% Hydromechanical wiht single-phase fluid flow finite element model.
%
% Each node has 3 degrees of freedom (dof):
% - 2 displacement components (ux,uy)
% - 1 pore pressure (p)
%
%% Author
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class definition
classdef Model_HM < Model_M     
    %% Constructor method
    methods
        function this = Model_HM()
            this = this@Model_M();
            this.ndof_nd = 3;       % Number of dofs per node
            this.physics = 'HM';    % Tag with the physics name
            disp("*** Physics: Hydromechanical with single-phase flow");
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
                        'fluid',this.mat.fluid);
                udofs = this.getElementDofs(el,[1,2]);
                pdofs = this.getElementDofs(el,3);
                elements(el) = RegularElement_HM(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,udofs,pdofs, ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                            this.isPlaneStress);
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
        % Plot the mesh with the boundary conditions
        function plotPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotPressureAlongSegment(Xi, Xf, npts, axisPlot);
        end

    end
        %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           ux        uy        Pl\n');
        end

    end
end