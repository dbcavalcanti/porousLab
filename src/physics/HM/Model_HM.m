%% Model_HM class
%
% Hydromechanical wiht single-phase fluid flow finite element model.
%
% Each node has 3 degrees of freedom (dof):
% - 2 displacement components (ux,uy)
% - 1 pore pressure (p)
%
%% Author
% Danilo Cavalcanti
%
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
            EFEMdraw.plotPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        %------------------------------------------------------------------
        % Update the result nodes coordinates of each element
        function updateResultVertices(this,configuration,factor)
            
            for el = 1:this.nelem
                
                % Initialize the vertices array
                vertices = this.element(el).type.result.vertices0;

                % Get the updated vertices:
                if strcmp(configuration,'Deformed')

                    % Update the nodal displacement vector associated to the
                    % element. This displacement can contain the enhancement
                    % degrees of freedom.
                    this.element(el).type.ue = this.U(this.element(el).type.gle); 

                    % Update the vertices based on the displacement vector
                    % associated to the element
                    for i = 1:length(this.element(el).type.result.faces)
                        X = vertices(i,:);
                        u = this.element(el).type.displacementField(X);
                        vertices(i,:) = X + factor*u';
                    end
                end
                this.element(el).type.result.setVertices(vertices);

            end
            
        end

        %------------------------------------------------------------------
        % Update the result nodes data of each element
        function updateResultVertexData(this,type)
            for el = 1:this.nelem
                % Update the nodal displacement vector associated to the
                % element. This displacement can contain the enhancement
                % degrees of freedom.
                this.element(el).type.ue = this.U(this.element(el).type.gle); 
                vertexData = zeros(length(this.element(el).type.result.faces),1);
                for i = 1:length(this.element(el).type.result.faces)
                    X = this.element(el).type.result.vertices(i,:);
                    if strcmp(type,'Model')
                        vertexData(i) = this.matID(el);
                    elseif strcmp(type,'Pressure')
                        p = this.element(el).type.pressureField(X);
                        vertexData(i) = p;
                    elseif strcmp(type,'Ux')
                        u = this.element(el).type.displacementField(X);
                        vertexData(i) = u(1);
                    elseif strcmp(type,'Uy')
                        u = this.element(el).type.displacementField(X);
                        vertexData(i) = u(2);
                    elseif strcmp(type,'E1')
                        s = this.element(el).type.strainField(X);
                        sp = this.element(el).type.principalStrain(s);
                        vertexData(i) = sp(1);
                    elseif strcmp(type,'PEMAG')
                        pe = this.element(el).type.plasticstrainMagnitude(X);
                        vertexData(i) = pe;
                    elseif strcmp(type,'Sx')
                        s = this.element(el).type.stressField(X);
                        vertexData(i) = s(1);
                    elseif strcmp(type,'Sy')
                        s = this.element(el).type.stressField(X);
                        vertexData(i) = s(2);
                    elseif strcmp(type,'Sxy')
                        s = this.element(el).type.stressField(X);
                        vertexData(i) = s(3);
                    elseif strcmp(type,'S1')
                        s = this.element(el).type.stressField(X);
                        sp = this.element(el).type.principalStress(s);
                        vertexData(i) = sp(1);
                    elseif strcmp(type,'S2')
                        s = this.element(el).type.stressField(X);
                        sp = this.element(el).type.principalStress(s);
                        vertexData(i) = sp(2);
                    elseif strcmp(type,'Sr')
                        s = this.element(el).type.stressField(X);
                        sp = this.element(el).type.stressCylindrical(s,X);
                        vertexData(i) = sp(1);
                    end
                end
                this.element(el).type.result.setVertexData(vertexData);
            end
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