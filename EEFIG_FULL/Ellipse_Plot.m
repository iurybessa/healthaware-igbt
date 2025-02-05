function Ellipse_Plot(A, C, varargin)
%
%  Ellipse_Plot(A,C,N) plots a 2D ellipse or a 3D ellipsoid 
%  represented in the "center" form:  
%               
%                   (x-C)' A (x-C) <= 1
%
%  Inputs: 
%  A: a 2x2 or 3x3 matrix.
%  C: a 2D or a 3D vector which represents the center of the ellipsoid.
%  PlotNo: A number that will be written on the plot
%
%  Originally Written by: Nima Moshtagh
%  Edited By: Masud Moshtaghi
%  m.moshtaghi@pgrad.unimelb.edu
%  University of Melbourne
%  Apr 20, 2009


%%%%%%%%%%%  start  %%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 50; % Default value for grid
PlotColor=rand(1,3);
% See if the user wants a different value for N.
%----------------------------------------------
if nargin >= 3
%   	N = varargin{1};
    PlotColor = varargin{1};
end
if nargin > 3
%   	N = varargin{1};
    PlotNo = varargin{2};
end

% check the dimension of the inputs: 2D or 3D
%--------------------------------------------
if length(C) == 3,
    Type = '3D';
elseif length(C) == 2,
    Type = '2D';
else
    display('Cannot plot an ellipse with more than 3 dimensions!!');
    return
end

% "singular value decomposition" to extract the orientation and the
% axes of the ellipsoid
[U D V] = svd(A);

if strcmp(Type, '2D'),
    % get the major and minor axes
    %------------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));

    theta = (0:1/N:2*pi+1/N);

    % Parametric equation of the ellipse
    %----------------------------------------
    state(1,:) = a*cos(theta); 
    state(2,:) = b*sin(theta);

    % Coordinate transform 
    %----------------------------------------
    X = V * state;
    X(1,:) = X(1,:) + C(1);
    X(2,:) = X(2,:) + C(2);
    
elseif strcmp(Type,'3D'),
    % generate the ellipsoid at (0,0,0)
    %----------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    c = 1/sqrt(D(3,3));
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);
    
    %  rotate and center the ellipsoid to the actual center point
    %------------------------------------------------------------
    XX = zeros(N+1,N+1);
    YY = zeros(N+1,N+1);
    ZZ = zeros(N+1,N+1);
    for k = 1:length(X),
        for j = 1:length(X),
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V * point;
            XX(k,j) = P(1)+C(1);
            YY(k,j) = P(2)+C(2);
            ZZ(k,j) = P(3)+C(3);
        end
    end
end


% Plot the ellipse
%----------------------------------------
if strcmp(Type,'2D'),
 
    plot(X(1,:),X(2,:),'-g','LineWidth',3,'Color',PlotColor);
    hold on;
    [D I]=max(X,[],2);
    if(nargin > 3)
        text(X(1,I(1,1)),X(2,I(1,1)),PlotNo,'FontSize',10,'FontWeight','bold');
    end
%     elseif(PlotNo<0)
%         text(X(1,I(1,1)),X(2,I(1,1)),[num2str(-PlotNo), '\sigma']);
%     end
    
    hold on;
    plot(C(1),C(2),'b*');
%     axis equal grid
else
    mesh(XX,YY,ZZ,'EdgeColor',PlotColor,'EdgeAlpha',0.3);
    axis equal
    hidden off
end

