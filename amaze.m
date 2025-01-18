function [B,BpX,BpY,sol,sx,sy] = amaze(n,branching)
    % A maze generator.
    %
    % [B,BpX,BpY,C,CpX,CpY] = amaze(n,branching)
    %
    % n, maze is n-by-n.
    % branching = 'first', 'middle', or 'last'.
    %    Experiment with all three.
    %
    % Defaults:
    % amaze(15,'middle')
    %
    % See also: amazing, graph, graph/shortestpath, delsq, numgrid,
    % https://blogs.mathworks.com/cleve/2019/03/04.
    
    % Copyright 2019 Cleve Moler
    % Copyright 2019 The MathWorks, Inc.
        
    if nargin < 1
        n = 15;
    end
    if nargin < 2
        branching = 'middle';
    end

    % Barriers, initially n^2 barriers.
    % Discrete Laplacian on a square grid.

    Db = delsq(numgrid('S',n+2));
    B = graph(logical(Db),'omitselfloops');

    % Cells, initially just nodes, no edges.
    m = n-1;  % m^2 cells
    Dc = delsq(numgrid('S',m+2));
    C = graph(logical(Dc),'omitselfloops');
        
    available = 1:m^2; % Nodes we haven't visited yet.
    branches = [];
    tree = zeros(0,2); % Depth first search.
    p = 1;

    while 1  % Break when available is empty
        available(p) = 0;
        if ~any(available)
            break
        end

        [~,~,ngh] = find(available(neighbors(C,p)));

        if ~isempty(ngh)
            idx = randi(length(ngh));  % Random choice.
            q = ngh(idx);              % Next cell.
            if length(ngh) > 1
                % Could have chosen another neighbor.
                branches(end+1) = p;
            end

            % Add a cell and remove a barrier.
            tree(end+1,:) = [p q];
            [i,j] = wall(p,q,m);
            B = rmedge(B,i,j);

            p = q;

        else
            
            for p = branches
                if all(available(neighbors(C,p)) == 0)
                    branches(branches==p) = [];
                end
            end

            % Take another branch.
            switch branching
                case 'first'
                    idx = 1;
                case 'last'
                    idx = length(branches);
                otherwise
                    idx = round(length(branches)/2);
            end

            p = branches(idx);
            branches(idx) = [];
        end
    end

    B = rmedge(B,1,n+1);
    B = rmedge(B,n^2-n,n^2);
    C = graph(tree(:,1),tree(:,2));
    [BpX,BpY,CpX,CpY] = makegrid;

    [solpath,len] = shortestpath(C,1,numnodes(C)); 
    sol = subgraph(C,solpath);
    sx = CpX(solpath);
    sy = CpY(solpath);


    % -------------------------------------------------------------

    function [i,j] = wall(p,q,m)
        % Wall [i,j] blocks path [p,q].
        switch q-p
            case -m  % west
                i = p+ceil(p/m)-1;
                j = i+1;
            case -1  % north
                i = p+ceil(p/m)-1;
                j = i+n;
            case 1  % south
                i = p+ceil(p/m);
                j = i+n;
            case m  % east
                i = p+ceil(p/m)-1+n;
                j = i+1;
        end
    end

    function [BpX,BpY,CpX,CpY] = makegrid       
        k = (0:n^2-1)';
        BpX = floor(k/n);
        BpY = mod(n-k-1,n);

        k = 1:m^2;
        CpX = floor((k-1)/m)+0.5;
        CpY = mod(n-k-1,m)+0.5;

    end
end