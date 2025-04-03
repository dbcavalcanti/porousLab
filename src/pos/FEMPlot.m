%% FEMDraw class
%
% This class implements methods to plot graphical results
% from the FEM analysis
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef FEMPlot < handle
    %% Public properties
    properties
        model = [];                     % Handle to an object of the model class
    end
    %% Constant properties
    properties (Constant, Access = private)
        DEFAULT_FONT_NAME = 'Times';    % Default font name for plots
        DEFAULT_FONT_SIZE = 16;         % Default font size for plots
        DEFAULT_COLOR = '-b';           % Default color and line style for plots
        DEFAULT_LINE_WIDTH = 1.5;       % Default line width for plots
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = FEMPlot(model)
            if (nargin > 0)
                this.model = model;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes the bounding box of the model with a tolerance.
        % Note:
        %   - The tolerance is calculated as 10% of the absolute difference 
        %     between the maximum and minimum x-coordinates of the model's nodes.
        %
        function bbox = getBoundingBox(this)
            minX = Inf; minY = Inf; maxX = -Inf; maxY = -Inf;
            for el = 1:this.model.nelem
                res = this.model.element(el).type.result;
                minX = min(minX,min(res.vertices(:,1)));
                minY = min(minY,min(res.vertices(:,2)));
                maxX = max(maxX,max(res.vertices(:,1)));
                maxY = max(maxY,max(res.vertices(:,2)));
            end
            tolx = 0.10*abs(max(this.model.NODE(:,1)) - min(this.model.NODE(:,1)));
            toly = 0.10*abs(max(this.model.NODE(:,1)) - min(this.model.NODE(:,1)));
            bbox(1) = minX - tolx;
            bbox(2) = minY - toly;
            bbox(3) = maxX + tolx;
            bbox(4) = maxY + toly;
        end

        %------------------------------------------------------------------
        function elements(this)
            
            % Initialize combined matrices
            allFaces      = [];      % Combined faces connectivity
            allVertices   = [];      % Combined vertices coordinates
            allVertexData = [];      % Combined vertex data
            vertexOffset  = 0;       % Offset for face indices due to combined vertices
            
            % Loop through all elements to build combined matrices
            for el = 1:this.model.nelem
                % Get the current element and its results
                element = this.model.element(el).type;
                res = element.result;
                
                % Adjust face indices for combined vertices
                faces = res.faces + vertexOffset;
                
                % Append to combined matrices
                allFaces      = [allFaces; faces];
                allVertices   = [allVertices; res.vertices];
                allVertexData = [allVertexData; res.vertexData];
                
                % Update vertex offset for the next element
                vertexOffset = vertexOffset + size(res.vertices, 1);
            end
            
            % Draw all patches at once
            patch('Faces', allFaces, ...
                'Vertices', allVertices, ...
                'FaceVertexCData', allVertexData, ...
                'FaceColor', 'interp', ...  % Use 'interp' for vertex-based colors
                'LineWidth', this.model.element(1).type.result.edgesThickness, ...
                'LineStyle', '-', ...
                'Marker', 'none');  % Disable markers for batch plotting
            
            % Set the colormap
            colormap(this.model.element(1).type.result.colormapType);
        end
        
        %------------------------------------------------------------------
        % Plots the continuum elements of the mesh
        function plotMesh(this)
            
            figure
            hold on, box on, grid on
            axis equal
            
            % Plot the continuum elements
            this.elements();
            
            % Get the bounding box
            bbox = this.getBoundingBox();
            xlim([bbox(1) bbox(3)]);
            ylim([bbox(2) bbox(4)]);
            
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            
        end
        
        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotGasPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            
            % Coordinates of the points where the pressure field is going
            % to be evaluated
            X = [linspace(Xi(1), Xf(1),npts);linspace(Xi(2), Xf(2),npts)]';
            
            % Initialize the pressure vector at these points
            P = zeros(size(X,1),1);
            S = zeros(size(X,1),1);
            
            % Calculate the pressure field in the points of the segment
            for i = 1:npts
                
                % Longitudinal coordinate of the point X(i) wrt to Xi
                S(i) = sqrt((X(i,1) - Xi(1))*(X(i,1) - Xi(1)) + (X(i,2) - Xi(2))*(X(i,2) - Xi(2)));
                
                % Find in each element this point is inside
                elem = findElementInMesh(this.model.NODE, this.model.ELEM, X(i,:));
                
                % Calculate the pressure field in the point X using the
                % shape function of the elem
                P(i) = this.model.element(elem).type.gasPressureField(X(i,:));
                
            end
            
            % Initialize, plot and configure the figure
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(P,S,'-b','LineWidth',1.5);
                xlabel('Gas pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,P,'-b','LineWidth',1.5);
                ylabel('Gas pressure (kPa)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            
        end
        
        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotCapillaryPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            
            % Coordinates of the points where the pressure field is going
            % to be evaluated
            X = [linspace(Xi(1), Xf(1),npts);linspace(Xi(2), Xf(2),npts)]';
            
            % Initialize the pressure vector at these points
            P = zeros(size(X,1),1);
            S = zeros(size(X,1),1);
            
            % Calculate the pressure field in the points of the segment
            for i = 1:npts
                
                % Longitudinal coordinate of the point X(i) wrt to Xi
                S(i) = sqrt((X(i,1) - Xi(1))*(X(i,1) - Xi(1)) + (X(i,2) - Xi(2))*(X(i,2) - Xi(2)));
                
                % Find in each element this point is inside
                elem = findElementInMesh(this.model.NODE, this.model.ELEM, X(i,:));
                
                % Calculate the pressure field in the point X using the
                % shape function of the elem
                P(i) = this.model.element(elem).type.capillaryPressureField(X(i,:));
                
            end
            
            % Initialize, plot and configure the figure
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(P,S,'-b','LineWidth',1.5);
                xlabel('Capillary pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,P,'-b','LineWidth',1.5);
                ylabel('Capillary pressure (kPa)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            
        end
        
        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            
            % Coordinates of the points where the pressure field is going
            % to be evaluated
            X = [linspace(Xi(1), Xf(1),npts);linspace(Xi(2), Xf(2),npts)]';
            
            % Initialize the pressure vector at these points
            P = zeros(size(X,1),1);
            S = zeros(size(X,1),1);
            
            % Calculate the pressure field in the points of the segment
            for i = 1:npts
                
                % Longitudinal coordinate of the point X(i) wrt to Xi
                S(i) = sqrt((X(i,1) - Xi(1))*(X(i,1) - Xi(1)) + (X(i,2) - Xi(2))*(X(i,2) - Xi(2)));
                
                % Find in each element this point is inside
                elem = findElementInMesh(this.model.NODE, this.model.ELEM, X(i,:));
                
                % Calculate the pressure field in the point X using the
                % shape function of the elem
                P(i) = this.model.element(elem).type.pressureField(X(i,:));
                
            end
            
            % Initialize, plot and configure the figure
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(P,S,'-b','LineWidth',1.5);
                xlabel('Pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,P,'-b','LineWidth',1.5);
                ylabel('Pressure (kPa)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            
        end
        
        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotDisplacementAlongSegment(this, dir, Xi, Xf, npts, axisPlot)
            
            % Coordinates of the points where the pressure field is going
            % to be evaluated
            X = [linspace(Xi(1), Xf(1),npts);linspace(Xi(2), Xf(2),npts)]';
            
            % Initialize the pressure vector at these points
            U = zeros(size(X,1),1);
            S = zeros(size(X,1),1);
            
            % Calculate the pressure field in the points of the segment
            for i = 1:npts
                
                % Longitudinal coordinate of the point X(i) wrt to Xi
                S(i) = sqrt((X(i,1) - Xi(1))*(X(i,1) - Xi(1)) + (X(i,2) - Xi(2))*(X(i,2) - Xi(2)));
                
                % Find in each element this point is inside
                elem = findElementInMesh(this.model.NODE, this.model.ELEM, X(i,:));
                
                % Calculate the pressure field in the point X using the
                % shape function of the elem
                u = this.model.element(elem).type.displacementField(X(i,:));
                U(i) = u(dir);
                
            end
            
            % Initialize, plot and configure the figure
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(U,S,'-b','LineWidth',1.5);
                if dir == 1
                    xlabel('Horizontal displacement (m)');
                elseif dir == 2
                    xlabel('Vertical displacement (m)');
                end
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,U,'-b','LineWidth',1.5);
                if dir == 1
                    ylabel('Horizontal displacement (m)');
                elseif dir == 2
                    ylabel('Vertical displacement (m)');
                end
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            
        end
        
        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotDisplacementAlongDiscontinuity(this, X0, axisPlot,resPath)
            
            mdl = this.model;
            
            % Initialize the stress vector at these points
            % wn = zeros(2*size(mdl.FRACT,1),1);
            % ws = zeros(2*size(mdl.FRACT,1),1);
            % S  = zeros(2*size(mdl.FRACT,1),1);
            
            count = 1;
            
            % Calculate the pressure field in the points of the segment
            for el = 1:mdl.nelem
                if sum(mdl.IDenr(el,:)) >= 1
                    for i = 1:mdl.element(el).type.fracture{1}.nIntPoints
                        
                        % Get the coordinates of the integration point
                        Xn = mdl.element(el).type.fracture{1}.intPoint(i).X;
                        X  = mdl.element(el).type.fracture{1}.shape.coordNaturalToCartesian(mdl.element(el).type.fracture{1}.node,Xn);
                        
                        % Longitudinal coordinate of the point X(i) wrt to X0
                        S(count) = norm(X-X0);
                        
                        % Get the stress vector
                        STRAIN = mdl.element(el).type.fracture{1}.intPoint(i).strain;
                        
                        % Save the stress component
                        ws(count) = STRAIN(1);
                        wn(count) = STRAIN(2);
                        
                        % Update counter
                        count = count + 1;
                        
                    end
                    
                end
                
            end
            
            [S,id] = sort(S);
            ws = ws(id);
            wn = wn(id);
            
            % Initialize, plot and configure the figure of the Shear stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(ws,S,'-b','LineWidth',1.5);
                xlabel('Tangential displacement (m)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,ws,'-b','LineWidth',1.5);
                ylabel('Tangential displacement (m)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
            if isempty(resPath)== false
                figName = fullfile(resPath,'faultTangentialDisplacement.fig');
                savefig(gcf,figName);
            end
            
            % Initialize, plot and configure the figure of the normal stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(wn,S,'-b','LineWidth',1.5);
                xlabel('Normal displacement (m)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,wn,'-b','LineWidth',1.5);
                ylabel('Normal displacement (m)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            if isempty(resPath)== false
                figName = fullfile(resPath,'faultNormalDisplacement.fig');
                savefig(gcf,figName);
            end
            
        end
        
        %------------------------------------------------------------------
        % Plot the cohesive stress field along a segment
        function plotCohesiveStressesAlongDiscontinuity(this, X0 ,axisPlot,resPath)
            
            mdl = this.model;
            
            % Initialize the stress vector at these points
            % tn    = zeros(2*size(mdl.FRACT,1),1);
            % ts    = zeros(2*size(mdl.FRACT,1),1);
            % S     = zeros(2*size(mdl.FRACT,1),1);
            % Saux  = zeros(2*size(mdl.FRACT,1),1);
            
            count = 1;
            
            % Calculate the pressure field in the points of the segment
            for el = 1:mdl.nelem
                if sum(mdl.IDenr(el,:)) >= 1
                    for i = 1:mdl.element(el).type.fracture{1}.nIntPoints
                        
                        % Get the coordinates of the integration point
                        Xn = mdl.element(el).type.fracture{1}.intPoint(i).X;
                        X  = mdl.element(el).type.fracture{1}.shape.coordNaturalToCartesian(mdl.element(el).type.fracture{1}.node,Xn);
                        
                        % Get the centroid of the fracture
                        Xref = mdl.element(el).type.fracture{1}.Xref;
                        
                        % Apply a pertubation to the coordinate X to avoid
                        % the duplicated points in S. It will "attrack" the
                        % point to the fracture centroid
                        % dX = 1.0e-5*(Xref-X)/norm(Xref-X);
                        
                        % Longitudinal coordinate of the point X(i) wrt to X0
                        S(count)    = norm(X-X0);
                        % Saux(count) = norm(X + dX -X0);
                        
                        % Get the stress vector
                        STRESS = mdl.element(el).type.fracture{1}.intPoint(i).stress;
                        
                        % Save the stress component
                        ts(count) = STRESS(1);
                        tn(count) = STRESS(2);
                        
                        % Update counter
                        count = count + 1;
                        
                    end
                    
                end
                
            end
            
            % Sort the vectors
            % [~,id] = sort(Saux);
            [~,id] = sort(S);
            S  = S(id);
            ts = ts(id);
            tn = tn(id);
            
            % Initialize, plot and configure the figure of the Shear stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(ts/1000,S,'-b','LineWidth',1.5);
                xlabel('Shear cohesive stress (MPa/m)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,ts/1000,'-b','LineWidth',1.5);
                ylabel('Shear cohesive stress (MPa/m)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            if isempty(resPath) == false
                figName = fullfile(resPath,'faultShearCohesiveStress.fig');
                savefig(gcf,figName);
            end
            
            % Initialize, plot and configure the figure of the normal stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(tn/1000,S,'+-b','LineWidth',1.5);
                xlabel('Normal cohesive stress (MPa/m)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,tn/1000,'+-b','LineWidth',1.5);
                ylabel('Normal cohesive stress (MPa/m)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontName',this.DEFAULT_FONT_NAME);
            set(gca,'FontSize',this.DEFAULT_FONT_SIZE);
            if isempty(resPath) == false
                figName = fullfile(resPath,'faultNormalCohesiveStress.fig');
                savefig(gcf,figName);
            end
            
        end
    end
end