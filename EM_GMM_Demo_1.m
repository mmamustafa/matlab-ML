% Expectation Maximization for Gaussian Mixture Models
% Reference:
%     http://www.ics.uci.edu/~smyth/courses/cs274/notes/EMnotes.pdf
%     http://www.cse.psu.edu/~rtc12/CSE586Spring2010/papers/prmlMixturesEM.pdf
% ----------------------------------------------------
%
% Variables:
%   <d>     space dimension
%   <N>     total number of data points
%   <n>     index of data point
%   <K>     total number of clusters
%   <k>     index of cluster
%   <C>     Cell array (1 x K): each cell is Clusters (data structure):
%           <u>     (d x 1)     mean
%           <S>     (d x d)     covariance
%           <p>     (1 x 1)     mixing coefficient


clc, clear all
toPlot = 1;     % works only for 1D, 2D, 3D

% DATA GENERATION
% ---------------
Dim = 2;                            % number of dimensions
NC = 4;                             % number of clusters
nPerCluster = 100;                  % number of points per cluster
rangeOfCentroidsPerDim = [-3 3];    % slection space for each dimension
WL = [-5 5 -5 5];

pts = [];                           % data points initialization
for i=1:NC
    % mean
    u = uniSampleND(repmat(rangeOfCentroidsPerDim,1,Dim));
    % Covariance
    while(1)
        t1 = -1+2*rand(Dim);    % random matrix with values between -1 and 1
        S = t1'*t1;
        cond(S)
        break
    end
    pts = [pts mvnrnd(u,S,nPerCluster)'];   % Generate from Multivariate Gaussians
    
    C_true{i}.u = u;
    C_true{i}.S = S;
end


% EXPECTATION MAXIMIZATION
% ------------------------
[C,gamma] = gmm_em(pts,NC,0,10);



% PLOTTING
% --------
if toPlot && Dim<4
    figure(1),    clf,
    set(gcf, 'color', [1 1 1])
    plotPoints(pts,'r*');
    hold on
    for k=1:NC
        h=error_ellipse(C{k}.S,C{k}.u,0.99);
        set(h,'Color','k')
    end
    hold off
    axis equal; axis(WL);
end


