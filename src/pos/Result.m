%% Result Class
% This class defines a _Result_ object to use the patch function from 
% Matlab to plot the model.
%
%% Methods
% This class provides the following methods:
%
% * *Result*: Constructor that initializes a Result object with the given vertices, faces, vertex data, and data label.
% * *setDataLabel*: Sets the label of the result being plotted.
% * *setVertices*: Updates the vertices of the patches.
% * *setFaces*: Updates the faces of the patches.
% * *setVertexData*: Updates the values associated with the vertices.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00: Initial version (December 2022).
% 
% Initially prepared for the course CIV 2801 - Fundamentos de Computação
% Gráfica, 2022, second term, Department of Civil Engineering, PUC-Rio.
%
%% Class Definition
classdef Result < handle

    %% Public attributes
    properties
        % Patch data
        vertices0        = [];          % Orginal vertices of the patches
        vertices        = [];           % Vertices of the patches
        faces           = [];           % Faces definitions
        vertexData      = [];           % Values associated to the vertices
        dataLabel       = '';           % Label of the result being plotted

        % Patch properties
        faceColor       = 'interp';     % Face color scheme for plotting
        markerType      = 'None';          % Type of marker 
        markerFaceColor = 'k';          % Color of the face of the marker
        edgesThickness  = 0.7;          % Thickness of the lines
        % colormapType    = 'jet(13)';        % Type of colormap scheme
        colormapType    = 'jet';        % Type of colormap scheme
    end

    %% Constructor method
    methods
        function this = Result(vertices,faces,vertexData,dataLabel)
            if nargin > 0
                this.vertices0 = vertices;
                this.vertices = vertices;
                this.faces = faces;
                this.vertexData = vertexData;
                this.dataLabel = dataLabel;
            end
        end        
    end

    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Set the label of the result being plotted
        function setDataLabel(this, dataLabel)
            this.dataLabel = dataLabel;
        end

        % -----------------------------------------------------------------
        % Set the vertices of the patches
        function setVertices(this, vertices)
            this.vertices = vertices;
        end

        % -----------------------------------------------------------------
        % Set the faces of the patches
        function setFaces(this, faces)
            this.faces = faces;
        end

        % -----------------------------------------------------------------
        % Set the values of vertices of the patches
        function setVertexData(this, vertexData)
            this.vertexData = vertexData;
        end

    end
end