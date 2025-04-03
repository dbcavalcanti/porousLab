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
        DEFAULT_FONT_NAME           = 'Times';        % Default font name for plots
        DEFAULT_FONT_SIZE           = 16;             % Default font size for plots
        DEFAULT_COLOR               = '-b';           % Default color and line style for plots
        DEFAULT_LINE_WIDTH          = 1.5;            % Default line width for plots
        DEFAULT_FACE_COLOR          = 'interp';       % Default face color scheme for plotting
        DEFAULT_MARKER_TYPE         = 'None';         % Default type of marker
        DEFAULT_MARKER_FACE_COLOR   = 'k';            % Default color of the face of the marker
        DEFAULT_EDGES_THICKNESS     = 0.7;            % Default thickness of the lines
        DEFAULT_COLORMAP_TYPE       = 'jet';          % Default colormap scheme
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
                'FaceColor', this.DEFAULT_FACE_COLOR, ...  % Use 'interp' for vertex-based colors
                'LineWidth', this.DEFAULT_EDGES_THICKNESS, ...
                'LineStyle', '-', ...
                'Marker', this.DEFAULT_MARKER_TYPE);  % Disable markers for batch plotting
            
            % Set the colormap
            colormap(this.DEFAULT_COLORMAP_TYPE);
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
                plot(P,S,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
                xlabel('Gas pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,P,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
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
                plot(P,S,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
                xlabel('Capillary pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,P,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
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
                plot(P,S,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
                xlabel('Pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,P,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
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
                plot(U,S,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
                if dir == 1
                    xlabel('Horizontal displacement (m)');
                elseif dir == 2
                    xlabel('Vertical displacement (m)');
                end
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,U,this.DEFAULT_COLOR,'LineWidth',this.DEFAULT_LINE_WIDTH);
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
    end
end