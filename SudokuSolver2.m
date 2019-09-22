% Sudoku Solver 9X9
% Solve by backtracking
%
% Mohamed Mustafa,  January 2016

function SudokuSolver2()
clc

% (1) Initial graph, each cell has array of all possible values.
%     Cell (0,0) is the lower left cell.
%     All cells are initialized with -1. (Cell can accept any value)
NC = 81;    % number of cells
graph = -ones(1,NC);

% (2) Draw the grid and accept initial known clues
fid = 1;
[graph,flag] = editGraph(graph,fid);

% load graph.mat
% graph = cell2mat(graph);

plotGrid(graph)

% (3) Check if graph is solvable...

% (4) Solve by backtracking
tic
[graph,status] = solver(graph,0);
toc

% (5) plot the result
plotGrid(graph)

return
% -------------------------------------------------------------------------
function plotGrid(graph,fid)

if nargin<2
    fid = 1;
end
figure(fid), clf

% (1) plot axis
for i=0:9
    % line width
    if mod(i,3)==0,     ms = 5;
    else                ms = 1;     end
    % Horizontal lines
    plot([0 9],[i i],'k-','linewidth',ms)
    if i==0,
        hold on,    set(gca,'YTick',[]);    set(gca,'XTick',[]);
    end
    % Vertical lines
    plot([i i],[0 9],'k-','linewidth',ms)
    set(gcf, 'color', [1 1 1])
    axis equal,axis([0 9 0 9])
end

% (2) add known values only if they have one VALID value (1-9)
for i=1:length(graph)
    if graph(i)>0
        [x,y] = ind2xy(i);
        text(x,y,num2str(graph(i)),'FontSize',30,'horizontalalignment','center')
    end
end
hold off
return
% -------------------------------------------------------------------------
function [graph,flag] = editGraph(graph,fid)
% flag:     0 --> known >= 17
%           1 --> knowns < 17 (less than minimum --> terminate)

if nargin<2
    fid = 1;
end

% Plot initial graph
plotGrid(graph,fid);
stop = 0;   % stopping condition
disp('*** Editting initial Sudoku Grid. ***')
disp('    >> Press Esc to finish.')
disp('    >> Press Space to clear cell.')
while ~stop
    % Get mouse-click position
    [x,y,key] = ginput(1);
    t1 = str2double(char(key));     % change key format
    % edit the graph
    if t1>0 && t1<10
        i = xy2ind(x,y);
        graph(i) = t1;
    end
    % clear curret cell
    if key==32
        i = xy2ind(x,y);
        graph(i) = -1;
    end
    % plot the updated graph
    plotGrid(graph,fid)
    % stopping condition
    if key==27
        stop = 1;
    end
end

% Check the number of known cells
c = 0;
for i=1:length(graph)
    if graph(i)>0
        c = c+1;
    end
end
flag = 0;
if c<17
    disp('*** Number of known clues is LESS than 17! ***')
    flag = 1;
end
return
% -------------------------------------------------------------------------
function [x,y] = ind2xy(i)
x = mod(i-1,9) + 0.5;
y = floor((i-1)/9) + 0.5;
return
% -------------------------------------------------------------------------
function i = xy2ind(x,y)
i = floor(y)*9 + floor(x) + 1;
return
% -------------------------------------------------------------------------
function out = getHorizontalSingles(graph,x,y)
x = setdiff(1:9,ceil(x))-0.5;
i = xy2ind(x,repmat(y,1,8));
out = [];
for j=1:8
    t1 = graph(i(j));
    if length(t1)==1 && t1(1)>0
        out = [out t1];
    end
end
out = unique(out);
return
% -------------------------------------------------------------------------
function out = getVerticalSingle(graph,x,y)
y = setdiff(1:9,ceil(y))-0.5;
i = xy2ind(repmat(x,1,8),y);
out = [];
for j=1:8
    t1 = graph(i(j));
    if length(t1)==1 && t1(1)>0
        out = [out t1];
    end
end
out = unique(out);
return
% -------------------------------------------------------------------------
function out = getSquareSingles(graph,x,y)
t1 = ceil(x/3);     x = (1:3)+(t1-1)*3;
t2 = ceil(y/3);     y = (1:3)+(t2-1)*3;
[x,y] = meshgrid(x,y);
x = x(:)-0.5;       y = y(:)-0.5;
i = xy2ind(x,y);
out = [];
for j=1:9
    t1 = graph(i(j));
    if length(t1)==1 && t1(1)>0
        out = [out t1];
    end
end
out = unique(out);
return
% -------------------------------------------------------------------------
function [graph,status] = solver(graph,toPlot)
% recursive solver.
if nargin<2
    toPlot = 0;
end

% (1) find next unassigned cell. If empty --> done
t0 = find(graph==-1);
if isempty(t0)
    status = 1;
    return
end
ind = t0(1);

% (2) find all possible/valid values
[x,y] = ind2xy(ind);
t1 = 1:sqrt(length(graph));
out = getHorizontalSingles(graph,x,y);
t1 = setdiff(t1,out);
out = getVerticalSingle(graph,x,y);
t1 = setdiff(t1,out);
out = getSquareSingles(graph,x,y);
t1 = setdiff(t1,out);

% (3) try each possible value (loop)
for possibleValue = t1
    % (4) update the graph
    graph(ind) = possibleValue;
    
    % plot for illustration
    if toPlot,      plotGrid(graph),    end
    
    % (5) call solver recurssively
    [graph1,status1] = solver(graph,toPlot);
    
    if status1
        graph = graph1;
        status = 1;
        return
    end
end

status = 0;
return



