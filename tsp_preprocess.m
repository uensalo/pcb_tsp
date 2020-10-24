%Author: Unsal Ozturk

clear
clc

COMPUTE_TC = 1;
COMPUTE_FW = 1; % NOTE: MAY TAKE LONG, ~3 to 4 MINUTES.

points = xlsread('data.xlsx', 'Points ');
rects = xlsread('data.xlsx', 'Blocks');
rects(:,3:4) = rects(:,1:2) + rects(:,3:4);

M = zeros(max(points));
n = size(M,1) * size(M,2); %number of nodes in the grid

for i = 1:length(rects)
    a = rects(i,1:2);
    b = rects(i,3:4);
    M(a(1):b(1), a(2):b(2)) = 1;
end

%COMPUTES PATH THROUGH SHORTEST MANHATTAN DISTANCES.
%DRILL MAY NOT MOVE TO ANOTHER HOLE IF SHORTEST MANHATTAN DISTANCE IS
%OBSTRUCTED...
%compute manhattan distances and paths
if(COMPUTE_TC)
fprintf('Computing distances with shortest Manhattan distance only movement for %d points and %d nodes.\n', length(points), n);
tic
dist = zeros(length(points) * (length(points) - 1) / 2, 1);
k = 1;
paths = cell(length(points) * (length(points) - 1) / 2, 1);
for i = 1:length(points) - 1
    for j = i + 1:length(points)
        path = shortest_taxicab_path(M, points(i,:), points(j,:));
        reachableij = ~isempty(path);
        paths{k} = path;
        if(reachableij)
            dist(k) = sum(abs(points(i,:) - points(j,:)));
        else
            dist(k) = Inf;
        end
        k = k + 1;
    end
end
fprintf('Done. ');
toc
end

%DRILL MAY MOVE TO ANY HOLE
%SHORTEST PATHS THAT ARE NOT THE SHORTEST MANHATTAN PATHS ALLOWED
%set up edge weight matrix based on adjacencies   
x = inf(n,n);
for i = 1:n
    %down
    x(i,i) = 0;
    if(i <= n - size(M,2))
        x(i, i + size(M,2)) = 1;
    end
    %up
    if(i > size(M,2))
        x(i, i - size(M,2)) = 1;
    end
    %right
    if(mod(i,size(M,2)) ~= 0)
        x(i, i + 1) = 1;
    end
    %left
    if(mod(i,size(M,2)) ~= 1)
        x(i, i - 1) = 1;
    end
end

for r = 1 : length(rects)
    aa = rects(r,1:2);
    bb = rects(r,3:4);
    for k = aa(1) : bb(1)
        for l = aa(2) : bb(2)
            i = coord2idx(size(M,2), k, l);
            if(i <= n - size(M,2))
                x(i, i + size(M,2)) = inf;
            end
            if(i > size(M,2))
                x(i, i - size(M,2)) = inf;
            end
            if(mod(i,size(M,2)) ~= 0)
                x(i, i + 1) = inf;
            end
            if(mod(i,size(M,2)) ~= 1)
                x(i, i - 1) = inf;
            end
        end
    end
end

if(COMPUTE_FW)
fprintf('Computing all pairs shortest paths for any drill movement for %d points and %d nodes.\nFloyd-Warshall for %d nodes may take some time.\n', length(points), n, n);
tic
[d_fw, next_fw] = floyd_warshall(x);
dist_fw = inf(length(points) * (length(points) - 1) / 2, 1);
paths_fw = cell(length(points) * (length(points) - 1) / 2, 1);
k = 1;
fprintf('Enumerating shortest paths...\n');
for i = 1:length(points) - 1
    for j = i + 1:length(points)
        pi = points(i,:);
        pj = points(j,:);
        idx_i = coord2idx(size(M,2),pi(1),pi(2));
        idx_j = coord2idx(size(M,2),pj(1),pj(2));
        path_fw = floyd_warshall_path(idx_i, idx_j, next_fw);
        path_converted = zeros(d_fw(idx_i,idx_j),2);
        for q = 1 : length(path_fw)
            [row_p, col_p] = idx2coord(path_fw(q), size(M,2));
            path_converted(q,:) = [row_p, col_p];
        end
        paths_fw{k} = path_converted;
        dist_fw(k) = d_fw(idx_i,idx_j);
        k = k + 1;
    end
end
fprintf('Done. ');
toc
end

fprintf('Writing paths and distances as tables...\n');
if(COMPUTE_FW)
table_fw_paths = cell2table(paths_fw(1:end,:));
writetable(table_fw_paths, 'table_fw_paths.xlsx');
writetable(table_fw_paths, 'table_fw_paths.csv');
end

if(COMPUTE_TC)
table_tc_paths = cell2table(paths(1:end,:));
writetable(table_tc_paths, 'table_tc_paths.xlsx');
writetable(table_tc_paths, 'table_tc_paths.csv');
end

if(COMPUTE_TC)
table_tc_dist = array2table(dist);
writetable(table_tc_dist, 'table_tc_dist.xlsx');
writetable(table_tc_dist, 'table_tc_dist.csv');
end

if(COMPUTE_FW)
table_fw_dist = array2table(dist_fw);
writetable(table_fw_dist, 'table_fw_dist.xlsx');
writetable(table_fw_dist, 'table_fw_dist.csv');
end

writetable(array2table(points), 'points.csv');
writetable(array2table(rects), 'rects.csv');

function path = shortest_taxicab_path(M, src, dst)
    dstx = dst(1);
    dsty = dst(2);
    con = 0;
    %initialize visited matrix
    visited = zeros(size(M));
    %initialize stack and push the source node
    s = src;
    prevx = -1;
    prevy = -1;
    parent = cell(size(M));
    while(~isempty(s))
        %pop from stack
        curx = s(1,1);
        cury = s(1,2);
        s(1,:) = []; 
        %if node is a block, don't even visit the node.
        if(M(curx,cury) == 1)
            continue;
        elseif(curx == dstx && cury == dsty)
            %destination found, return true
            con = 1;
            parent{curx,cury} = [prevx, prevy];
            break;
        elseif(visited(curx,cury) ~= 1)
            %mark as visited
            visited(curx,cury) = 1;
            %record parent node
            parent{curx,cury} = [prevx, prevy]; 
            prevx = curx;
            prevy = cury;
            %compute neighbours to push to stack
            dirx = 0;
            if(curx > dstx)
                dirx = -1;
            elseif(curx < dstx)
                dirx = 1;
            end
            diry = 0;
            if(cury > dsty)
                diry = -1;
            elseif(cury < dsty)
                diry = 1;
            end         
            %push neighbours down the stack
            s = [curx + dirx, cury; s];
            s = [curx, cury + diry; s];
        end
    end   
    %reconstruct path from source to destination using parent pointers
    if(con)
        curx = dstx;
        cury = dsty;
        taxicab_dist = sum(abs(src - dst));
        k = taxicab_dist + 1;
        path = zeros(k,2);
        path(k,:)= [dstx, dsty];
        k = k - 1;
        while(k > 0)        
            par =  parent{curx,cury};
            curx = par(1);
            cury = par(2);
            path(k,:) = [curx, cury];
            k = k - 1;
        end
    else
        path = [];
    end
end

function [dist, next] = floyd_warshall(x)
    n = length(x);
    dist = x;
    next = -ones(n);
    for i = 1 : n
        for j = 1 : n
            next(i,j) = j;
        end
    end
    for i = 1 : n
        dist(i,i) = 0;
        next(i,i) = i;
    end
    for k = 1 : n
        if(mod(k,ceil(n/10)) == 0)
            fprintf('\tFloyd-Warshall: %3.1f%% complete.\n', k / n * 100);
        end
        for i = 1 : n
            for j = 1 :n
                if(dist(i,j) > dist(i,k) + dist(k,j))
                    dist(i,j) = dist(i,k) + dist(k,j);
                    next(i,j) = next(i,k);
                end
            end
        end
    end
end

function path = floyd_warshall_path(src, dst, next)
    if(next(src,dst) == -1)
        path = [];
    else
        path = src;
        cur = src;
        while(cur ~= dst)
            cur = next(cur,dst);
            path = [path, cur];
        end
    end
end

function [row, col] = idx2coord(idx, c)
    col = mod(idx,c);
    if(int32(col) == 0)
        col = c;
    end
    row = ceil(idx / c);
end

function idx = coord2idx(cols, r, c)
    idx = (r - 1) * cols + c;
end
