function varargout = fmpmPotentialPDE2d(xdomain, ydomain, varargin)
%fmpmPotentialPDE2d computes the influence potential field of the
% force and moment partitioning method [1] using the MATLAB Partial
% Differential Equation Toolbox.
%
% Two-Dimensional version developed based on code used in Zhu et al. [2].
%
% The boundary conditions on the non-reference geometry are
% (n dot grad(phi)) = 0 in multi-body configurations, which implies the
% non-penetration property of the influence field.
% -------------------------------------------------------------------------
%
% Syntax:
%
% phi = fmpmPotentialPDE2d(xdomain, ydomain, geo, 'loadtype')
% phi = fmpmPotentialPDE2d(xdomain, ydomain, geo1, geo2, ..., geoN, 'loadtype')
% phi = fmpmPotentialPDE2d(__, refgeoN)
% phi = fmpmPotentialPDE2d(__, refgeoN, innerTrSize, outerTrSize)
% [phi, X, Y] = fmpmPotentialPDE2d(__)
% [__, bodygeo] = fmpmPotentialPDE2d(__)
% [__, bodygeo, x_surf, y_surf] = fmpmPotentialPDE2d(__)
% [__, bodygeo, x_surf, y_surf, phi_surf] = fmpmPotentialPDE2d(__)
% [__, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf] = fmpmPotentialPDE2d(__)
% [__, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf, norm_surf] = fmpmPotentialPDE2d(__)
% -------------------------------------------------------------------------
%
% Inputs:
% 
% xdomain: x-coordiante vector of the computational domain (outter boundary).
%
% ydomain: y-coordiante vector of the computational domain (outter boundary).
%
% geo1, geo2, ..., geoN: N-by-2 vectors give the x and y coordinates of the
%                        body geometries. (inner boundary). The first
%                        geometry geo1 is the reference geometry for which
%                        the phi field will be computed by default. The
%                        referece geometry/geometries can be specified by
%                        the refgeoN argument otherwise.
%
% 'loadtype': a char variable that specifies the type of load of which the
%             influence field will be computed. 'loadtype' must be 'lift',
%             'pitch', or 'drag'. When 'pitch' is used, a 1-by-2 array,
%             refaxis, that specifies the x and y coordinates of the pitch
%             moment reference axis should follow the 'pitch' argument. If
%             refaxis is empty [], the default refaxis is the geometric
%             center of the reference geometry.
%
% refgeoN: a scalar or a row vector that specifies the reference
%          geometry/geometries for the phi field. The values of the
%          elements refer to the order of geo1, geo2, ..., geoN inputs. The
%          default reference geometry is geo1 if refgeoN is not specified
%          or refgeoN = [];
%
% innerTrSize: a scalar that specifies the size of the triangular mesh on
%              the inner boudaries (bodies). The default value is 0.5*ds
%              where ds is the minimum unit length of the computational
%              domain if innerTrSize is not specified or innerTrSize = [];
%
% outerTrSize: a scalar that specifies the size of the triangular mesh on
%              the outer boudaries (computational domain). The default
%              value is 10*ds where ds is the minimum unit length of the
%              computational domain if outerTrSize is not specified or
%              outerTrSize = [];
%
% -------------------------------------------------------------------------
%
% Outputs:
%
% phi: influence potential field.
%
% X: x-coordinate mesh generated based on xdomain and ydomain.
%
% Y: y-coordinate mesh generated based on xdomain and ydomain.
%
% bodygeo: polyshape object defined by the geo input.
%
% x_surf: x-coordinate vector of the data points on the reference body
%         surface.
%
% y_surf: y-coordinate vector of the data points on the reference body
%         surface.
%
% phi_surf: phi vector of the data points on the reference body surface.
%
% dphidx_surf: vector of gradients of phi in the x-direction on the data
%              points on the reference body surface.
%
% dphidy_surf: vector of gradients of phi in the y-direction on the data
%              points on the reference body surface.
%
% normxy_surf: an N-by-2 array that contains the x- (first column) and y-
%              (second column) componets of the surface normal vectors of
%              the reference geometries. The dimension and the order of
%              normxy_surf corespond to x_surf and y_surf.
% -------------------------------------------------------------------------
%
% Examples:
%
% 1. geo = [[-0.5, 0.5, 0]', [-0.5, -0.5, 0.5]'];
%    [phi, X, Y] = fmpmPotentialPDE2d(-2:0.1:2, -1:0.1:1, geo, 'lift');
%
%    Computes the lift influence potential field on a triangle defiend by
%    points (-0.5, -0.5), (0.5, -0.5), and (0, 0.5) in a computational
%    domain defiend by x-coordinates -2:0.1:2 and y-coordinates -1:0.1:1.
%    Output phi returns the influence field; X and Y return the mesh grids.
%
% 2. geo1 = [[-0.5, 0.5, 0]', [-0.5, -0.5, 0.5]'];
%    geo2 = [[0, 1.2, 1.2, 0]', [-0.2, -0.2, 0.1, 0.1]'];
%    [phi, ~, ~, bodyshape] = ...
%    fmpmPotentialPDE2d(-2:0.1:2, -1:0.1:1, geo1, geo2, 'pitch', [0.5, 0], 2, 0.05, 0.1);
%
%    Computes the pitch moment influence potential field on geo2
%    (as refgeoN = 2) in a multi-body geometry by geo1 and geo2. The size
%    of the triangular mesh on the geometry is 0.05 on the inner boundaries
%    and 0.1 on the outer boudaries. The pitch moment reference point is
%    given by refaxis = [0.5, 0]. Output bodyshape returns the combined
%    geometry as a polyshape object.
%
% 3. geo1 = [[-0.5, 0.5, 0]', [-0.5, -0.5, 0.5]'];
%    geo2 = [[0, 1.2, 1.2, 0]', [-0.2, -0.2, 0.1, 0.1]'];
%    geo3 = [[0.1, 1.3, 1.3, 0.1]', [-0.3, -0.3, 0, 0]'];
%    phi = fmpmPotentialPDE2d(-2:0.1:2, -1:0.1:1, geo1, geo2, geo3, 'pitch', [], [1, 3]);
%
%    Computes the pitch moment influence potential field on the combined
%    body of geo1 and geo3 (as refgeoN = [1, 3]) referenced to the
%    geometric center (as refaxis = []) of the combined geometry.
%
% 4. geo1 = [[-0.5, 0.5, 0]', [-0.5, -0.5, 0.5]'];
%    geo2 = [[0, 1.2, 1.2, 0]', [-0.2, -0.2, 0.1, 0.1]'];
%    [phi, ~, ~, bodyshape, x_surf, y_surf, phi_surf, ~, ~, normxy_surf] = ...
%    fmpmPotentialPDE2d(-2:0.1:2, -1:0.1:1, geo1, geo2, 'lift');
%
%    Computes the lift influence potential field of geo1 (default setting
%    as refgeoN is not specified) and returns the polygon object of the
%    conbined geometry that contains geo1 and geo2. x_surf and y_surf
%    returns the coordinates of the surface data points on the reference
%    geometry geo1. phi_surf returns the corresponding phi values and
%    normxy_surf returns the corresponding surface normal vectors of these
%    points.
% -------------------------------------------------------------------------
%
% version 2.0.1
%   - Fixed a bug where using a single output argument causes a formatting
%     error.
% Xiaowei He and Yuanhang Zhu
% 7/20/2023
% -------------------------------------------------------------------------
% version 2.0.0
%   - Extended the influence field to multi-body configurations.
%   - Added an output argument for surface normal vectors on the reference
%     geometry.
%   - Optimized the organization of the input and output arguments
%     accordingly.
%   - Updates in the headlines and examples.
% Xiaowei He and Yuanhang Zhu
% 7/17/2023
% -------------------------------------------------------------------------
% version 1.2.1
%   - Optimized the outer boundary detection to prevent potential false
%     boundary conditions.
%   - Changed the default minimum triangular size to 0.5*innerTrSize to
%     prevent potential insufficient element size.
% Xiaowei He and Yuanhang Zhu
% 6/19/2023
% -------------------------------------------------------------------------
% version 1.1.0
%   - Added output options including phi values on the body surface and the
%     corresponding x- and y- gradients. See e.g.3 and e.g.4.
% Xiaowei He and Yuanhang Zhu
% 6/15/2023
% -------------------------------------------------------------------------
% version 1.0.0
% Xiaowei He and Yuanhang Zhu
% 05/30/2023
% -------------------------------------------------------------------------
%
% [1]. Menon, K. & Mittal, R. (2021). Quantitative analysis of the
% kinematics and induced aerodynamic loading of individual vortices in
% vortex-dominated flows: A computation and data-driven approach. Journal
% of Computational Physics 443, 110515.
%
% [2]. Zhu et al. (2023). Force moment partitioning and scaling analysis of
% vortices shed by a 2D pitching wing in quiescent fluid. arXiv:2301.13373
% [physics.flu-dyn].
% -------------------------------------------------------------------------

    nargoutchk(1, 10)
    
    % check output arguments
    flag_surf = false;
    switch nargout
        case 1
            outmsg = 'phi';
        case 3
            outmsg = '[phi, X, Y]';
        case 4
            outmsg = '[phi, X, Y, bodygeo]';
        case 6
            outmsg = '[phi, X, Y, bodygeo, x_surf, y_surf]';
            flag_surf = true;
        case 7
            outmsg = '[phi, X, Y, bodygeo, x_surf, y_surf, phi_surf]';
            flag_surf = true;
        case 9
            outmsg = '[phi, X, Y, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf]';
            flag_surf = true;
        case 10
            outmsg = '[phi, X, Y, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf, normxy_surf]';
            flag_surf = true;
        otherwise
            error('errorType:NumOutput', ...
                  ['Output variables must be in the form of' ...
                   '\n\tphi = fmpmPotentialPDE2d(__),' ...
                   '\n\t[phi, X, Y] = fmpmPotentialPDE2d(__) or,' ...
                   '\n\t[phi, X, Y, bodygeo] = fmpmPotentialPDE2d(__) or,' ...
                   '\n\t[phi, X, Y, bodygeo, x_surf, y_surf] = fmpmPotentialPDE2d(__) or,' ...
                   '\n\t[phi, X, Y, bodygeo, x_surf, y_surf, phi_surf] = fmpmPotentialPDE2d(__) or,' ...
                   '\n\t[phi, X, Y, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf] = fmpmPotentialPDE2d(__) or,' ...
                   '\n\t[phi, X, Y, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf, normxy_surf] = fmpmPotentialPDE2d(__).'])
    end
    
    % check input arguments
    if nargin < 4
        error('errorType:NumInput', ...
              [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''loadtype''):' ...
               '\n\tnot sufficient input arguments, the x- and y-corrdinates of the computational domain, at least one body geometry, and the ''loadtype'' must be specified.'])
    end

    % x and y domain
    if ~(isnumeric(xdomain) && isvector(xdomain) && isnumeric(ydomain) && isvector(ydomain))
        error('errorType:domain', ...
              [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''loadtype''):' ...
               '\n\txdomain and ydomain must be 1-D numeric vectors.'])
    end
    
    % body geometries
    iloadtype = 1;
    while ~ischar(varargin{iloadtype})
        iloadtype = iloadtype + 1;
    end
    geo = varargin(1:iloadtype-1);
    for i = 1 : length(geo)
        if ~(isnumeric(geo{i}) && size(geo{i}, 2) == 2)
            error('errorType:geo', ...
                  [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo1, geo2, ..., geoN, ''loadtype''):' ...
                   '\n\tgeo1, geo2, ..., geoN must be n-by-2 numeric matrices that contain the x and y coordiantes of the body geometries.'])
        end
    end
    
    % loadtype
    loadtype = varargin{iloadtype};
    if ~(strcmp(loadtype, 'lift') || strcmp(loadtype, 'pitch') || strcmp(loadtype, 'drag'))
        error('errorType:loadtype', ...
              [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''loadtype''):' ...
               '\n\t''loadtype'' must be ''lift'', ''pitch'', or ''drag''.'])
    end
    iTrSize = iloadtype+1;
    if strcmp(loadtype, 'pitch')
        iTrSize = iloadtype+2;
        if length(varargin) < iloadtype+1
            error('errorType:refaxis', ...
                  [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''pitch'', refaxis):' ...
                   '\n\trefaxis must be specified by an 1-by-2 numeric vector or an empty input [] when computing ''pitch'' field.'])
        else
            refaxis = varargin{iloadtype+1};
        end
        if ~(isnumeric(refaxis) && isrow(refaxis) && length(refaxis) == 2) && ~isempty(refaxis)
            error('errorType:refaxis', ...
                  [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''pitch'', refaxis):' ...
                   '\n\trefaxis must be specified by an 1-by-2 numeric vector or an empty input [] when computing ''pitch'' field.'])
        end
    end

    % tranglular sizes and refgeoN input
    innerTrSize = [];
    outerTrSize = [];
    refgeoN = 1;
    if length(varargin) >= iTrSize
        restIn = varargin(iTrSize:end);
        switch length(restIn)
            case 1
                refgeoN = restIn{1};
            case 2
                refgeoN = restIn{1};
                innerTrSize = restIn{2};
            case 3
                refgeoN = restIn{1};
                innerTrSize = restIn{2};
                outerTrSize = restIn{3};
            otherwise
                error('errorType:NumInput', ...
                      [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo1, geo2, ..., geoN, ''loadtype'', refgeoN, innerTrSize, outerTrSize):' ...
                       '\n\ttoo many input arguments.'])
        end
    end

    % triangular mesh parameters
    % unit lengths of x and y
    dx = abs(xdomain(2) - xdomain(1));
    dy = abs(ydomain(2) - ydomain(1));
    ds = min([dx, dy]);
    % default inner and outer triangle size
    if isempty(innerTrSize)
        innerTrSize = 0.5*ds;
    end
    if isempty(outerTrSize)
        outerTrSize = 10*ds;
    end
    if ~(isscalar(innerTrSize) && isnumeric(innerTrSize))
        error('errorType:TrSize', ...
              [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''loadtype'', refgeoN, innerTrSize, outerTrSize):' ...
               '\n\tinnerTrSize must be a numeric scalar or an empty input [].'])
    end
    if ~(isscalar(outerTrSize) && isnumeric(outerTrSize))
        error('errorType:TrSize', ...
              [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''loadtype'', refgeoN, innerTrSize, outerTrSize):' ...
               '\n\touterTrSize must be a numeric scalar or an empty input [].'])
    end

    % index numbers of reference geometries
    if isempty(refgeoN)
        refgeoN = 1;
    end
    if ~(isnumeric(refgeoN) && isrow(refgeoN))
        error('errorType:refgeoN', ...
              [outmsg ' = fmpmPotentialPDE2d(xdomain, ydomain, geo, ''loadtype'', refgeoN):' ...
               '\n\trefgeoN must be a numeric scalar, a row vector, or an empty input [].'])
    end

    % cartesian mesh of the domain
    [X, Y] = meshgrid(xdomain, ydomain);
    
    % geometry of the outer boundarys
    outer = polyshape([xdomain(1), xdomain(end), xdomain(end), xdomain(1)], ...
                      [ydomain(1), ydomain(1), ydomain(end), ydomain(end)]);
    
    % geometry of the inner boundarys
    NBody = length(geo);
    body = cell(1, NBody);
    bodyscaled = cell(1, NBody);
    inner = polyshape(geo{1});
    innerscaled = polyshape(geo{1});
    for i = 1 : NBody
        body{i} = polyshape(geo{i});
        [xc, yc] = centroid(body{i});
        bodyscaled{i} = scale(body{i}, 1.01, [xc, yc]);
        inner = union(inner, body{i});
        innerscaled = union(innerscaled, bodyscaled{i});
    end
    refgeo = body{refgeoN(1)};
    refgeoscaled = bodyscaled{refgeoN(1)};
    for i = 1 : length(refgeoN)
        refgeo = union(refgeo, body{refgeoN(i)});
        refgeoscaled = union(refgeoscaled, bodyscaled{refgeoN(i)});
    end
    if strcmp(loadtype, 'pitch') && isempty(refaxis)
        [xcref, ycref] = centroid(refgeo);
        refaxis = [xcref, ycref];
    end
    
    % computational domain
    domain = xor(outer, inner);
    
    % create pde model
    PDEmodel = createpde;
    % triangular mesh within the domain
    tr = triangulation(domain);
    tnodes = tr.Points';
    telements = tr.ConnectivityList';
    % create the mesh in the model
    geometryFromMesh(PDEmodel, tnodes, telements);
    hfig = figure('Visible', 'off');
    hEdge = pdegplot(PDEmodel,'EdgeLabels','on');
    hold on
    edgeX = hEdge.Parent.Children(1).VertexData(1, :)';
    edgeY = hEdge.Parent.Children(1).VertexData(2, :)';
    inInner = inpolygon(edgeX', edgeY', innerscaled.Vertices(:, 1), innerscaled.Vertices(:, 2));
    inRef = inpolygon(edgeX', edgeY', refgeoscaled.Vertices(:, 1), refgeoscaled.Vertices(:, 2));
    edgeInner = find(inInner);
    edgeOuter = setdiff(1:PDEmodel.Geometry.NumEdges, edgeInner);
    edgeRef = find(inRef);
    close(hfig)
    % adaptive triangular mesh
    generateMesh(PDEmodel, 'Hedge', {edgeOuter, outerTrSize, edgeInner, innerTrSize}, ...
                 'Hmax', 5*outerTrSize, 'Hmin', 0.5*innerTrSize, 'GeometricOrder', 'quadratic');
    % surface norm initialization
    surfnorm = [];
    surfx = [];
    surfy = [];
    % boundary conditions
    applyBoundaryCondition(PDEmodel, 'neumann', 'Edge', edgeOuter, 'g', 0, 'q', 0);
    applyBoundaryCondition(PDEmodel, 'neumann', 'Edge', edgeRef, 'g', @bcRefgeo, 'q', 0);
    if length(geo) > 1
        applyBoundaryCondition(PDEmodel, 'neumann', 'Edge', setdiff(edgeInner, edgeRef), 'g', 0, 'q', 0);
    end
    % solve pde
    specifyCoefficients(PDEmodel, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);
    result = solvepde(PDEmodel);
    phi = griddata(result.Mesh.Nodes(1, :), result.Mesh.Nodes(2, :), result.NodalSolution, X, Y) ...
                   - mean(result.NodalSolution);

    % phi on the body surface
    if flag_surf
        % obtain free boundaries
        [~, e, ~] = result.Mesh.meshToPet();
        ik1 = e(1, :);
        ik2 = e(2, :);
        x_surf = [result.Mesh.Nodes(1, ik1)'; result.Mesh.Nodes(1, ik2)'];
        y_surf = [result.Mesh.Nodes(2, ik1)'; result.Mesh.Nodes(2, ik2)'];
        phi_surf = [result.NodalSolution(ik1); result.NodalSolution(ik2)];
        dphidx_surf = [result.XGradients(ik1); result.XGradients(ik2)];
        dphidy_surf = [result.YGradients(ik1); result.YGradients(ik2)];
        % obtain the data on the refgeo surface
        inRefSurf = inpolygon(x_surf, y_surf, refgeoscaled.Vertices(:, 1), refgeoscaled.Vertices(:, 2));
        x_surf = x_surf(inRefSurf);
        y_surf = y_surf(inRefSurf);
        phi_surf = phi_surf(inRefSurf) - mean(result.NodalSolution);
        dphidx_surf = dphidx_surf(inRefSurf);
        dphidy_surf = dphidy_surf(inRefSurf);
        normx = griddata(surfx, surfy, surfnorm(:, 1), x_surf, y_surf, 'nearest');
        normy = griddata(surfx, surfy, surfnorm(:, 2), x_surf, y_surf, 'nearest');
        normxy_surf = [normx, normy];
    end

    % output
    % [phi, X, Y, bodygeo, x_surf, y_surf, phi_surf, dphidx_surf, dphidy_surf, normxy_surf]
    switch nargout
        case 1
            varargout{1} = phi;
        case 3
            varargout{1} = phi;
            varargout{2} = X;
            varargout{3} = Y;
        case 4
            varargout{1} = phi;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = inner;
        case 6
            varargout{1} = phi;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = inner;
            varargout{5} = x_surf;
            varargout{6} = y_surf;
        case 7
            varargout{1} = phi;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = inner;
            varargout{5} = x_surf;
            varargout{6} = y_surf;
            varargout{7} = phi_surf;
        case 9
            varargout{1} = phi;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = inner;
            varargout{5} = x_surf;
            varargout{6} = y_surf;
            varargout{7} = phi_surf;
            varargout{8} = dphidx_surf;
            varargout{9} = dphidy_surf;
        case 10
            varargout{1} = phi;
            varargout{2} = X;
            varargout{3} = Y;
            varargout{4} = inner;
            varargout{5} = x_surf;
            varargout{6} = y_surf;
            varargout{7} = phi_surf;
            varargout{8} = dphidx_surf;
            varargout{9} = dphidy_surf;
            varargout{10} = normxy_surf;
    end
    
    % nonconstant bc
    function bcMatrix = bcRefgeo(location, state)
        switch loadtype
            case 'lift'
                bcMatrix = location.ny;
            case 'pitch'
                bcMatrix = dot(cross([location.x-refaxis(1), location.y-refaxis(2), 0], ...
                                     [location.nx, location.ny, 0]), [0 0 1], 2);
            case 'drag'
                bcMatrix = location.nx;
            otherwise
                error('''loadtype'' must be ''lift'', ''pitch'', or ''drag''.')
        end
        surfnorm = [surfnorm; [location.nx, location.ny]];
        surfx = [surfx; location.x];
        surfy = [surfy; location.y];
    end
end