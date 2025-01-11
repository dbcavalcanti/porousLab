%% EFEMDraw class
%
% This class implements methods to plot graphical results
% from the EFEM analysis
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef EFEMDraw < handle
    %% Public properties
    properties
        model = []; % handle to an object of the EFEMmodel class
    end
    
    %% Class (constant) properties
    properties (Constant)
        supsize_fac = 0.035; % factor for support size
        loadsize_fac = 0.3; % factor for maximum load size in relation to
                             % half vertical window size
        minloadsize_fac = 0.012 % factor for minimum load size
        arrowsize_fac = 0.025; % factor for load arrow size
        loadstep_fac = 0.05; % factor for load step
        dimlineshift_fac = 0.85; % factor for down shift of dimension lines
                                 % with respect to half Y size
        dimlinetick_fac = 0.15; % factor for dimension line tick sizeEFEM
                                % with respect to half Y size
        maxdisplsize_fac = 0.60; % factor for maximum transversal size in
                                 % relation to half vertical window size
        inflectpt_fac = 0.005; % factor for size of inflection point disk
        rotmeter_fac = 0.85; % factor for rotation meter size in
                             % relation to half vertical window size
        picktol_fac = 0.01; % factor for picking a point
        minmemblen_fac = 0.05; % factor for minimum member length
        ValidSupInsertion = 1; % status for valid support insertion
        BeamLineNotFound = 2; % status for beam line not found for support insertion
        SupInsertionNotValid = 3; % status for not valid position for support insertion
        SupDelMinNumSup = 1; % status for minimum number of internal supports
        ValidSupDeletion = 2; % status for valid support deletion
        SupDelNotFound = 3; % status for support not found for deletion
        MembLoadFound = 1; % status for pick member load found
        MembLoadNotFound = 2 % status for pick member not found
        SupMoveFound = 1; % status for support found for moving
        SupMoveNotFound = 2; % status for support not found for moving
    end
    
    %% Private properties
    properties (Access = private)
        pickmember = 0; % current member load picked
        picksup = 0; % current support picked
        orig_suppos = 0; % original moving support position
        hnd_draft = []; % handle to draft graphics object being displayed
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = EFEMDraw(model)
            if (nargin > 0)
                draw.model = model;
            end
        end
    end
    
    %% Class (static) auxiliary functions
    methods (Static)
        %------------------------------------------------------------------
        % Plots a square with defined center coordinates, side length and
        % color.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  l: side length
        %  c: color (RGB vector)
        function square(cnv,x,y,l,c)
            X = [x - l/2, x + l/2, x + l/2, x - l/2];
            Y = [y - l/2, y - l/2, y + l/2, y + l/2];
            fill(cnv, X, Y, c);
        end
        
        %------------------------------------------------------------------
        % Plots a draft version of a square with defined center coordinates,
        % side length and color.
        % It draws a hollow square and sets its parent property as the
        % given handle to the group of graphics objects.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  hnd: handle to group of graphics objects
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  l: side length
        %  c: color (RGB vector)
        function draftSquare(cnv,hnd,x,y,l,c)
            X = [x - l/2, x + l/2, x + l/2, x - l/2, x - l/2];
            Y = [y - l/2, y - l/2, y + l/2, y + l/2, y - l/2];
            plot(cnv, X, Y, 'color', c, 'Parent', hnd);
        end
        
        %------------------------------------------------------------------
        % Plots a triangle with defined top coordinates, height, base,
        % orientation, and color.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: triangle height
        %  b: triangle base
        %  ang: angle (in radian) between the axis of symmetry and the
        %       horizontal direction (counterclockwise) - 0 rad when
        %       triangle is pointing left
        %  c: color (RGB vector)
        function triangle(x,y,h,b,ang,c)
            cx = cos(ang);
            cy = sin(ang);

            X = [x, x + h * cx + b/2 * cy, x + h * cx - b/2 * cy];
            Y = [y, y + h * cy - b/2 * cx, y + h * cy + b/2 * cx];
            fill(X, Y, c);
        end
        
        %------------------------------------------------------------------
        % Plots a draft version of a triangle with defined top coordinates,
        % height, base, orientation, and color.
        % It draws a hollow triangle and sets its parent property as the
        % given handle to the group of graphics objects.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  hnd: handle to group of graphics objects
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: triangle height
        %  b: triangle base
        %  ang: angle (in radian) between the axis of symmetry and the
        %       horizontal direction (counterclockwise) - 0 rad when
        %       triangle is pointing left
        %  c: color (RGB vector)
        function draftTriangle(hnd,x,y,h,b,ang,c)
            cx = cos(ang);
            cy = sin(ang);

            X = [x, x + h * cx + b/2 * cy, x + h * cx - b/2 * cy, x];
            Y = [y, y + h * cy - b/2 * cx, y + h * cy + b/2 * cx, y];
            plot(X, Y, 'color', c, 'Parent', hnd);
        end
        
        %------------------------------------------------------------------
        % Plots a circle with defined center coordinates, radius and color.
        % This method is used to draw hinges on 2D models.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  r: circle radius
        %  c: color (RGB vector)
        function circle(cnv,x,y,r,c)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            plot(cnv, xcirc, ycirc, 'color', c);
        end
        
        %------------------------------------------------------------------
        % Plots a circle disk with defined center coordinates, radius and
        % color. The circle is filled with the given color
        % This method is used to inflection on 2D models.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  r: circle radius
        %  c: color (RGB vector)
        function disk(cnv,x,y,r,c)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            fill(cnv, xcirc, ycirc, c);
        end
        
        %------------------------------------------------------------------
        % Plots an arrow with defined beggining coordinates, length,
        % arrowhead height, arrowhead base, orientation, and color.
        % This method is used to draw load symbols on 2D models.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: beggining coordinate on the X axis
        %  y: beggining coordinate on the Y axis
        %  l: arrow length
        %  h: arrowhead height
        %  b: arrowhead base
        %  ang: pointing direction (angle in radian with the horizontal
        %       direction - counterclockwise) - 0 rad when pointing left
        %  c: color (RGB vector)
        function arrow2D(cnv,x,y,l,h,b,ang,c)
            cx = cos(ang);
            cy = sin(ang);

            X = [x, x + l * cx];
            Y = [y, y + l * cy];
            line(cnv, X, Y, 'Color', c,'LineWidth',2);
            EFEMDraw.triangle(x, y, h, b, ang, c);
        end

        %------------------------------------------------------------------
        % Snap a value to the closest step value.
        function snap_val = snapToStepValue(val,step)
            fp = val / step;   % "fraction" part
            ip = floor(fp);    % integer part
            fp = fp - ip;
            if fp > 0.5
                snap_val = (ip + 1.0) * step;
            elseif fp < -0.5
                snap_val = (ip - 1.0) * step;
            else
                snap_val = ip * step;
            end
        end
    end

    %% Protect methods
    methods (Access = public) % Access from methods in subclasses
        
        function max_load = getMaxLoad(draw)
            max_load = max(draw.model.F);
        end

        function elements(draw,plotMatId)
            nmaterials = size(unique(draw.model.matID),1);
            colors = parula(nmaterials);
            for el = 1:draw.model.nelem
                res = draw.model.element(el).type.result;
                if plotMatId
                    materialID = draw.model.matID(el);
                    fieldFaceColor = colors(materialID,:);
                else
                    fieldFaceColor = res.faceColor;
                end
                patch('Faces',res.faces,...
                    'Vertices',res.vertices,...
                    'FaceVertexCData',res.vertexData,...
                    'FaceColor',fieldFaceColor,...
                    'LineWidth',res.edgesThickness,...
                    'LineStyle','-',...
                    'Marker',res.markerType,...
                    'MarkerFaceColor',res.markerFaceColor); 
                colormap(res.colormapType);
                if isa(draw.model.element(el).type,'EnrichedElement')
                    res = draw.model.element(el).type.fracture{1}.result;
                    patch('Faces',res.faces,...
                        'Vertices',res.vertices,...
                        'FaceColor','white',...
                        'EdgeColor','k',...
                        'EdgeAlpha',0.5,...
                        'LineStyle','--',...
                        'LineWidth',res.edgesThickness,...
                        'Marker',res.markerType,...
                        'MarkerFaceColor',res.markerFaceColor); 
                end

            end
        end

        function fractures(draw)
            for el = 1:size(draw.model.FRACT,1)
                x = [draw.model.NODE_D(draw.model.FRACT(el,:),1)];
                y = [draw.model.NODE_D(draw.model.FRACT(el,:),2)];
                p = plot(x,y,'k--','MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1);
                % p.Color(4) = 0.3;
            end
        end
    end
    
    %% Public methods
    methods        
        %------------------------------------------------------------------
        % Returns the bounding box (x and y limits) of a continuous beam
        % model. The returned box has xmin = 0, xmax = totalLen,
        % ymin = -totalLen * 0.05, ymax = +totalLen * 0.05, in which
        % totalLen is the length of the entire beam model.
        % The y limits are fictitious. They are equal in module and 
        % equal to a small percentage of the total length to force 
        % the adjustment of the box in the y direction, keeping y = 0
        % in the center of the canvas.
        function bbox = EFEMBoundBox(draw)
            minX = Inf; minY = Inf; maxX = -Inf; maxY = -Inf;
            for el = 1:draw.model.nelem
                res = draw.model.element(el).type.result;
                minX = min(minX,min(res.vertices(:,1)));
                minY = min(minY,min(res.vertices(:,2)));
                maxX = max(maxX,max(res.vertices(:,1)));
                maxY = max(maxY,max(res.vertices(:,2)));
            end
            tolx = 0.10*abs(max(draw.model.NODE(:,1)) - min(draw.model.NODE(:,1)));
            toly = 0.10*abs(max(draw.model.NODE(:,1)) - min(draw.model.NODE(:,1)));
            bbox(1) = minX - tolx;
            bbox(2) = minY - toly;
            bbox(3) = maxX + tolx;
            bbox(4) = maxY + toly;
        end
        
        %------------------------------------------------------------------
        % Draws a continuous beam model with applied loads.
        % Input:
        % - cnv: graphics context (owning canvas)
        function mesh(draw,plotMatId)
            if nargin == 1
                plotMatId = false;
            end

            figure
            hold on, box on, grid on
            axis equal
            
            % Draw continuous mesh
            draw.elements(plotMatId);

            % Draw fractures
            % draw.fractures();

            % Get the bounding box
            bbox = draw.EFEMBoundBox();
            xlim([bbox(1) bbox(3)]);
            ylim([bbox(2) bbox(4)]);

            set(gca,'FontName','Times','FontSize',16)
            
        end

        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotPressureAlongSegment(draw, Xi, Xf, npts,axisPlot)

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
                elem = findElementInMesh(draw.model.NODE, draw.model.ELEM, X(i,:));

                % Calculate the pressure field in the point X using the
                % shape function of the elem
                P(i) = draw.model.element(elem).type.pressureField(X(i,:));

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
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');

        end

        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotDisplacementAlongSegment(draw, dir, Xi, Xf, npts, axisPlot)

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
                elem = findElementInMesh(draw.model.NODE, draw.model.ELEM, X(i,:));

                % Calculate the pressure field in the point X using the
                % shape function of the elem
                u = draw.model.element(elem).type.displacementField(X(i,:));
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
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');

        end

        %------------------------------------------------------------------
        % Plot the pressure field along a segment
        function plotDisplacementAlongDiscontinuity(draw, X0, axisPlot,resPath)

            mdl = draw.model;

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
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
            if isempty(resPath)== false
                figName = fullfile(resPath,'faultNormalDisplacement.fig');
                savefig(gcf,figName);
            end

        end

        %------------------------------------------------------------------
        % Plot the cohesive stress field along a segment
        function plotCohesiveStressesAlongDiscontinuity(draw, X0 ,axisPlot,resPath)

            mdl = draw.model;

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
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
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
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
            if isempty(resPath) == false
                figName = fullfile(resPath,'faultNormalCohesiveStress.fig');
                savefig(gcf,figName);
            end

        end

        %------------------------------------------------------------------
        % Plot the cohesive stress field along a segment
        function plotCohesiveStressesIncrementAlongDiscontinuity(draw, X0 ,axisPlot,resPath)

            mdl = draw.model;

            % Initialize the stress vector at these points
            tn = zeros(2*size(mdl.FRACT,1),1);
            ts = zeros(2*size(mdl.FRACT,1),1);
            S  = zeros(2*size(mdl.FRACT,1),1);
            Saux  = zeros(2*size(mdl.FRACT,1),1);

            count = 1;

            % Calculate the pressure field in the points of the segment
            for el = 1:mdl.nelem
                if sum(mdl.IDenr(el,:)) >= 1
                    for i = 1:mdl.element(el).type.fracture{1}.nIntPoints

                        % Get the coordinates of the integration point
                        Xn = mdl.element(el).type.fracture{1}.intPoint(i).X;
                        X  = mdl.element(el).type.fracture{1}.shape.coordNaturalToCartesian(mdl.element(el).type.fracture{1}.node,Xn);

                        SGeo = mdl.element(el).type.fracture{1}.geostaticStressAtPoint(X, mdl.element(el).type, mdl.K0, mdl.yTop);

                        % Get the centroid of the fracture
                        Xref = mdl.element(el).type.fracture{1}.Xref;

                        % Apply a pertubation to the coordinate X to avoid
                        % the duplicated points in S. It will "attrack" the
                        % point to the fracture centroid
                        dX = 1.0e-5*(Xref-X)/norm(Xref-X);

                        % Longitudinal coordinate of the point X(i) wrt to X0
                        S(count)    = norm(X-X0);
                        Saux(count) = norm(X + dX -X0);
        
                        % Get the stress vector
                        STRESS = mdl.element(el).type.fracture{1}.intPoint(i).stress;

                        % Save the stress component
                        ts(count) = STRESS(1) - SGeo(1);
                        tn(count) = STRESS(2) - SGeo(2);
                        
                        % Update counter
                        count = count + 1;

                    end

                end

            end

            % Sort the vectors
            [~,id] = sort(Saux);
            S  = S(id);
            ts = ts(id);
            tn = tn(id);

            % Initialize, plot and configure the figure of the Shear stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(ts,S,'-b','LineWidth',1.5);
                xlabel('Shear cohesive stress increment (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,ts,'-b','LineWidth',1.5);
                ylabel('Shear cohesive stress increment (kPa)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
            figName = fullfile(resPath,'faultShearCohesiveStressIncrement.fig');
            savefig(gcf,figName);

            % Initialize, plot and configure the figure of the normal stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(tn,S,'-b','LineWidth',1.5);
                xlabel('Normal cohesive stress increment (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,tn,'-b','LineWidth',1.5);
                ylabel('Normal cohesive stress increment (kPa)');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
            figName = fullfile(resPath,'faultNormalCohesiveStressIncrement.fig');
            savefig(gcf,figName);

        end

        %------------------------------------------------------------------
        % Plot the cohesive stress field along a segment
        function plotSlipTendencyAlongDiscontinuity(draw, X0 ,axisPlot)

            mdl = draw.model;

            % Initialize the stress vector at these points
            ST = zeros(2*size(mdl.FRACT,1),1);
            S  = zeros(2*size(mdl.FRACT,1),1);

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
                        ST(count) = mdl.element(el).type.fracture{1}.intPoint(i).statevar(end);
                        
                        % Update counter
                        count = count + 1;

                    end

                end

            end

            [S,id] = sort(S);
            ST = ST(id);

            % Initialize, plot and configure the figure of the Shear stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(ST,S,'o-b','LineWidth',1.5);
                xlabel('Slip tendency');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,ST,'o-b','LineWidth',1.5);
                ylabel('Slip tendency');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');

        end

        %------------------------------------------------------------------
        % Plot the cohesive stress field along a segment
        function plotConformityErrorAlongDiscontinuity(draw, err, X0 ,axisPlot)

            mdl = draw.model;

            % Number of fracture nodes
            nfracNodes = size(mdl.NODE_D,1);

            % Initialize the path vector
            S  = zeros(nfracNodes,1);

            % Longitudinal coordinate of the point X(i) wrt to X0
            for i = 1:nfracNodes
                X  = mdl.NODE_D(i,:);
                S(i) = norm(X-X0);
            end

            [S,id] = sort(S);
            err = err(id);

            % Initialize, plot and configure the figure of the Shear stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(err,S,'o-b','LineWidth',1.5);
                xlabel('Conformity error');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,err,'o-b','LineWidth',1.5);
                ylabel('Conformity error');
                xlabel('Longitudinal distance (m)');
            end
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');

        end

        %------------------------------------------------------------------
        % Plot the cohesive stress field along a segment
        function plotPressureAlongDiscontinuity(draw, X0 ,axisPlot)

            mdl = draw.model;

            % Initialize the stress vector at these points
            P     = zeros(2*size(mdl.FRACT,1),1);
            Ptop  = zeros(2*size(mdl.FRACT,1),1);
            Pbot  = zeros(2*size(mdl.FRACT,1),1);
            S     = zeros(2*size(mdl.FRACT,1),1);
            Saux  = zeros(2*size(mdl.FRACT,1),1);

            count = 1;

            % Calculate the pressure field in the points of the segment
            for el = 1:mdl.nelem
                if mdl.element(el).type.isEnriched
                    Ue = mdl.U(mdl.element(el).type.gle);
                    % Number of pressure dofs associated with the continuum
                    nglpc = mdl.element(el).type.nglptot - mdl.element(el).type.nglpenrPf;
                    if mdl.element(el).type.staticCondensation.DisplJump
                        pc = mdl.element(el).type.ue(1+mdl.element(el).type.nglu:this.nglu+nglpc);
                    else
                        pc = mdl.element(el).type.ue(1+mdl.element(el).type.nglutot:mdl.element(el).type.nglutot+nglpc);
                    end 
                    pf = mdl.element(el).type.getDiscontinuityMidPlanePressure(pc);
                    pc = Ue(1+mdl.element(el).type.nglutot:mdl.element(el).type.nglutot+nglpc);
                    for i = 1:mdl.element(el).type.fracture{1}.nIntPoints

                        % Get the coordinates of the integration point
                        Xn = mdl.element(el).type.fracture{1}.intPoint(i).X;
                        X  = mdl.element(el).type.fracture{1}.shape.coordNaturalToCartesian(mdl.element(el).type.fracture{1}.node,Xn);
                        
                        % Get the centroid of the fracture
                        Xref = mdl.element(el).type.fracture{1}.Xref;

                        % Apply a pertubation to the coordinate X to avoid
                        % the duplicated points in S. It will "attrack" the
                        % point to the fracture centroid
                        dX = 1.0e-5*(Xref-X)/norm(Xref-X);

                        % Longitudinal coordinate of the point X(i) wrt to X0
                        S(count)    = norm(X-X0);
                        Saux(count) = norm(X + dX -X0);
        
                        % Get the mid-plane pressure
                        P(count) = pf(i);

                        % Natural coordinates associated with the continuum element
                        % of this point
                        Xnc = mdl.element(el).type.shape.coordCartesianToNatural(mdl.element(el).type.node,X);

                        % Shape function matrix of the continuum
                        N = mdl.element(el).type.shape.shapeFncMtrx(Xnc);

                        % Get the top plane pressure
                        Ntop = mdl.element(el).type.topEnhancedShapeFncMtrx(N,Xnc,1);
                        Ptop(count) = [N , Ntop]*pc;

                        % Get the bottom plane pressure
                        Nbot = mdl.element(el).type.bottomEnhancedShapeFncMtrx(N,Xnc,1);
                        Pbot(count) = [N , Nbot]*pc;
                        
                        % Update counter
                        count = count + 1;

                    end

                end

            end
            
            % Sort the vectors
            [~,id] = sort(Saux);
            S = S(id);
            P = P(id);
            Ptop = Ptop(id);
            Pbot = Pbot(id);

            % Initialize, plot and configure the figure of the Shear stress
            figure
            hold on, box on, grid on
            if strcmp(axisPlot,'y')
                plot(Ptop,S,'--k','LineWidth',1.5,'DisplayName','(+)-plane');
                plot(P,S,'-b','LineWidth',1.5,'DisplayName','Mid-plane');      
                plot(Pbot,S,'-.k','LineWidth',1.5,'DisplayName','(-)-plane');
                xlabel('Discontinuity internal pressure (kPa)');
                ylabel('Longitudinal distance (m)');
            elseif strcmp(axisPlot,'x')
                plot(S,Ptop,'--k','LineWidth',1.5,'DisplayName','(+)-plane');
                plot(S,P,'-b','LineWidth',1.5,'DisplayName','Mid-plane');       
                plot(S,Pbot,'-.k','LineWidth',1.5,'DisplayName','(-)-plane');
                ylabel('Discontinuity internal pressure (kPa)');
                xlabel('Longitudinal distance (m)');
            end
            legend;
            set(gca,'FontSize',16);
            set(gca,'FontName','Times');
            

        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a EFEMmodel object.
        function draw = clean(draw)
            draw.model = [];
        end
    end
end