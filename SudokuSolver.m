% Sudoku Solver 9X9
% Solve by searching/elimination
%
% Mohamed Mustafa,  December 2015

function SudokuSolver()
clc

% (1) Initial graph, each cell has array of all possible values.
%     Cell (0,0) is the lower left cell.
%     Empty cell --> No solution --> terminate
%     All cells are initialized with -1. (Cell can accept any value)
NC = 81;    % number of cells
graph = cell(1,NC);
for i=1:NC
    graph{i} = -1;
end

% % (2) Draw the grid and accept initial known clues
% fid = 1;
% [graph,flag] = editGraph(graph,fid);
% if flag
%     disp('*** Program terminated. ***')
%     return
% end
%
% save('graph.mat','graph'),return

load graph.mat
plotGrid(graph)

graph'

return

% (3) Eliminate all invalid value from each cell.
for i=1:NC
    t1 = graph{i};  % vector of possible values in that cell
    if isempty(t1)
        disp('*** No solution found! ***')
        return
    end
    % initialize empty cells
    if t1(1)<0
        t1 = 1:9;
    end
    % remove all invalid values
    if length(t1)>1
        % Remove from t1 vector
        [x,y] = ind2xy(i);
        out = getHorizontalSingles(graph,x,y);
        t1 = setdiff(t1,out);
        out = getVerticalSingle(graph,x,y);
        t1 = setdiff(t1,out);
        out = getSquareSingles(graph,x,y);
        t1 = setdiff(t1,out);
        % Update Graph
        graph{i} = t1;
        % if result is singleton, do something...
        
    end
end
plotGrid(graph)

% (4) Solve by backtracking...


return
% -------------------------------------------------------------------------
function plotGrid(graph,allCells,fid)
% allCells:     0  -->  only the singltons
%               1  -->  all vectors in the graph

if nargin<3
    fid = 1;
    if nargin<2
        allCells = 1;   % print all cells
    end
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
    t1 = length(graph{i});
    if t1==1
        if graph{i}>0
            [x,y] = ind2xy(i);
            text(x,y,num2str(graph{i}),'FontSize',30,'horizontalalignment','center')
        end
    end
    if t1>1 && allCells
        [x,y] = ind2xy(i);
        x0 = x-0.5+0.25;    y0 = y-0.5+0.25;
        x = x0;             y = y0;
        for j=1:t1
            text(x,y,num2str(graph{i}(j)),'FontSize',10,'horizontalalignment','center')
            x = x + 0.25;
            if x>(x0-0.25+0.9)
                x = x0;     y = y+0.25;
            end
        end
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
plotGrid(graph,1,fid);
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
        graph{i} = t1;
    end
    % clear curret cell
    if key==32
        i = xy2ind(x,y);
        graph{i} = -1;
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
    if graph{i}>0
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
    t1 = graph{i(j)};
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
    t1 = graph{i(j)};
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
    t1 = graph{i(j)};
    if length(t1)==1 && t1(1)>0
        out = [out t1];
    end
end
out = unique(out);
return
% -------------------------------------------------------------------------
function [graph,status] = solver(graph)
% recursive solver.

% (1) find next unassigned cell. If empty --> done



% (2) find all possible/valid values

% (3) try each possible value (loop)

return



