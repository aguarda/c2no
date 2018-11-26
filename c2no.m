function [varargout] = c2no(infile, outpcname, outfolder, ppc, varargin)
% c2no  Perform a Constant-size, Compact and Non-Overlapping clustering of a Point Cloud.
%
%       This function receives a Point Cloud, and outputs a set of clusters
%       with the same number of points, specified by ppc. If the total
%       number of points is not multiple of ppc, the user can specify if
%       only 1 cluster should be incomplete with the remaining points, or
%       if all points should be equally distributed. These incomplete
%       clusters can be padded by repeating points.
%       The resulting clusters are compact and do not overlap each other.
%       A detailed description of the algorithm is available at:
%
%   usage:  
%       c2no(infile, outpcname, outfolder, ppc)
%       c2no(infile, outpcname, outfolder, ppc, cnn)
%       c2no(infile, outpcname, outfolder, ppc, cnn, nb)
%       c2no(infile, outpcname, outfolder, ppc, cnn, nb, mic)
%       c2no(infile, outpcname, outfolder, ppc, cnn, nb, mic, padding)
%       [progress, clust_disp, clust_bb, sum_dists, iter_adj_rest] = c2no(infile, outpcname, outfolder, ppc, cnn, nb, mic, padding, debug)
%
%   input arguments:
%       infile - Input filename of the point cloud to cluster (PLY or PCD)
%       outpcname - Output clusters basename
%       outfolder - Output clusters directory path
%       ppc - Target number of points per cluster
%
%   optional input arguments:
%       cnn - Number of nearest neighboring clusters eligible for transfers
%           and refinement (Default: 10)
%       nb - Number of block divisions on each direction for Cluster
%           Centroid Initialization step, which divides the point cloud
%           into nb*nb*nb blocks (Default: 8)
%       mic - Selects the approach for dealing with target cluster sizes
%           when the total number of points in the point cloud is not
%           multiple of ppc; TRUE uses the Minimally Incomplete Clusters
%           case, where points are equally distributed through all
%           clusters; FALSE uses the Single Incomplete Cluster case, where
%           all but one cluster will have ppc points, and one cluster will
%           have the remaining points (Default: TRUE)
%       padding - Selects if padding is applied to incomplete clusters,
%           when the total number of points in the point cloud is not
%           multiple of ppc (Default: TRUE)
%       debug - This option allows tracking the progress of the algorithm
%           in the different stages, computing cluster dispersion metrics
%           at each iteration, and saving the clusters at the end of each
%           phase (Default: FALSE)
%
%   output arguments:
%       progress - 2D Matrix with the number of points in each cluster for
%           each iteration of the Iterative Clustering Phase
%       clust_disp - 3D Matrix with the cluster dispersion (measured as the
%           variance of cluster coordinates) for each of the 3 dimensions,
%           for each cluster, at each iteration of the Iterative Clustering
%           Phase (only if debug=TRUE)
%       clust_bb - 3D Matrix with the cluster dispersion (measured as the
%           range of cluster coordinates) for each of the 3 dimensions,
%           for each cluster, at each iteration of the Iterative Clustering
%           Phase (only if debug=TRUE)
%       sum_dists - 1D Matrix with the sum of distances between points and
%           respective centroids, averaged for all clusters, taken at each
%           iteration of the Clustering Refinement Phase (only if
%           debug=TRUE)
%       iter_adj_rest - Indicates at which iteration of the Iterative
%           Clustering Phase the adjacency restriction was lifted, or the
%           value 0 if it was not (only if debug=TRUE)
%
%   Author: André Guarda
%   email: andre.guarda@lx.it.pt
%   Release date: 23/11/2018

    % Check number of input parameters
    numvarargs = length(varargin);
    if numvarargs > 5
        error('c2no:TooManyInputs', 'Too many input arguments');
    end

    % Set default input arguments
    max_iter = 100;
    optargs = {10 8 true true false};

    % Overwrite default values with input
    optargs(1:numvarargs) = varargin;
    [num_nn, num_blocks, mic, padding, debug] = optargs{:};

    % Read the original point cloud file
    pc = pcread(infile);
    coordinates = unique(pc.Location,'rows'); % Discard repeated points
    % Determine number of clusters
    num_clust = ceil(length(coordinates)/ppc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        Initialization Phase                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %                  Cluster Centroids Initialization                  %
    
    centroids = [];
    box_tick = max([pc.XLimits(2) - pc.XLimits(1),...
                    pc.YLimits(2) - pc.YLimits(1),...
                    pc.ZLimits(2) - pc.ZLimits(1)])/num_blocks;
    oct = zeros(num_blocks,num_blocks,num_blocks);
    % Divide point cloud into blocks and count the points in each block
    for x = 0:num_blocks-1
        for y = 0:num_blocks-1
            for z = 0:num_blocks-1
                temp_box = coordinates(coordinates(:,1) >= (pc.XLimits(1) + x*box_tick)...
                                    & (coordinates(:,1) < (pc.XLimits(1) + (x+1)*box_tick) | x == num_blocks-1)...
                                    & coordinates(:,2) >= (pc.YLimits(1) + y*box_tick)...
                                    & (coordinates(:,2) < (pc.YLimits(1) + (y+1)*box_tick) | y == num_blocks-1)...
                                    & coordinates(:,3) >= (pc.ZLimits(1) + z*box_tick)...
                                    & (coordinates(:,3) < (pc.ZLimits(1) + (z+1)*box_tick) | z == num_blocks-1),:);
                
                oct(x+1,y+1,z+1) = length(temp_box);
            end
        end
    end
    % Attribute centroids to blocks with more points than PPC
    N = floor(oct./ppc);
    % Determine number of centroids left to attribute
    missing = num_clust - sum(sum(sum(N)));
    % Remove points (multiple of PPC) that already resulted in attributed centroids
    rem = mod(oct,ppc);
    % Sort blocks by number of remaining points (< PPC), and attribute 1
    % centroid to the ones with most points, until all centroids are
    % attributed
    [~,rem_idx] = sort(reshape(rem,1,[]),'descend');
    N(rem_idx(1:missing)) = N(rem_idx(1:missing)) + 1;
    % Now, compute centroid coordinates in each block
    for x = 0:num_blocks-1
        for y = 0:num_blocks-1
            for z = 0:num_blocks-1
                % Check the number of centroids attributed to this block
                temp_N = N(x+1,y+1,z+1);
                if temp_N == 0
                    continue;
                end
                % Get points of current block
                temp_box = coordinates(coordinates(:,1) >= (pc.XLimits(1) + x*box_tick)...
                                    & (coordinates(:,1) < (pc.XLimits(1) + (x+1)*box_tick) | x == num_blocks-1)...
                                    & coordinates(:,2) >= (pc.YLimits(1) + y*box_tick)...
                                    & (coordinates(:,2) < (pc.YLimits(1) + (y+1)*box_tick) | y == num_blocks-1)...
                                    & coordinates(:,3) >= (pc.ZLimits(1) + z*box_tick)...
                                    & (coordinates(:,3) < (pc.ZLimits(1) + (z+1)*box_tick) | z == num_blocks-1),:);
                
                new_octants = {temp_box};
                % Recursively divide the block in half, until all
                % sub-blocks have less than PPC points
                done_flag = false;
                while(~done_flag)
                    done_flag = true;
                    old_octants = new_octants;
                    new_octants = {};
                    for i=1:length(old_octants)
                        seg = old_octants{i};
                        % Decide if the sub-block is divided or not
                        if length(seg) > ppc
                            % Determine bounding box and where to cut
                            bbox = [max(seg); min(seg)];
                            [wid, idx] = max(bbox(1,:) - bbox(2,:));
                            cut_line = wid/2 + bbox(2,idx);
                            % Divide sub-block in half
                            new_octants{end + 1} = seg(seg(:,idx) < cut_line,:);
                            new_octants{end + 1} = seg(seg(:,idx) >= cut_line,:);

                            done_flag = false;
                        else
                            new_octants{end + 1} = old_octants{i};
                        end
                    end
                end
                % Sort all sub-blocks by number of points
                [~,oct_idx] = sort(cellfun('length', new_octants), 'descend');
                % Compute centroids for sub-blocks with most points
                for i = 1:temp_N
                    centroids = [centroids; single(mean(new_octants{oct_idx(i)}))];
                end
            end
        end
    end
    
    %                     Initial Cluster Assignment                     %
    
    [~,old_assignment] = pdist2(centroids, coordinates,'euclidean','Smallest',1);
    old_assignment = old_assignment';
    % Check if any cluster is empty
    bad_cent = setdiff(1:num_clust, unique(old_assignment));
    if ~isempty(bad_cent)
        % For each empty cluster, assign the closest point to it instead
        for i = 1:length(bad_cent)
            [~,c] = pdist2(coordinates, centroids(bad_cent(i),:),'euclidean','Smallest',1);
            old_assignment(c) = bad_cent(i);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     Iterative Clustering Phase                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check in total number of points of the point cloud is multiple of PPC
    % Set max and min cluster sizes for the MIC or SIC cases
    extra = mod(length(coordinates), ppc);
    idx_smaller_cluster = [];
    if extra == 0
        num_smaller = 0;
        max_size_cluster = ppc;
        min_size_cluster = ppc;
    elseif mic
        num_smaller = num_clust - mod(length(coordinates), num_clust);
        max_size_cluster = ceil(length(coordinates)/num_clust);
        min_size_cluster = floor(length(coordinates)/num_clust);
    else
        num_smaller = 1;
        max_size_cluster = ppc;
        min_size_cluster = extra
    end
    % Track number of points in each cluster at each iteration
    progress = zeros(num_clust,max_iter + 1);
    progress(:,1) = hist(old_assignment,num_clust)';
    % Track progress of algorithm for debugging purposes
    if debug
        % Iteration in which adjacency restriction is lifted
        iter_adj_rest = 0;
        % Variance and bounding box of each cluster at each iteration
        clust_disp = zeros(num_clust,3,max_iter + 2);
        clust_bb = zeros(num_clust,3,max_iter + 2);
        for j=1:num_clust
            cluster = coordinates(old_assignment == j, :);
            clust_disp(j,:,1) = var(cluster);
            clust_bb(j,:,1) = range(cluster);
        end
    end
    % Compute threshold for adjacency restriction %
    % Sample 1% of the points
    rem_idx = randsample(1:length(coordinates), round(0.01*length(coordinates)));
    current_points_samp = coordinates(rem_idx, :);
    % Compute distance to closest point for each
    [min_neigh_dist,~] = pdist2(coordinates,current_points_samp,'euclidean','Smallest',2);
    min_neigh_dist(1,:) = [];
    % Compute median and multiply by diagonal of unit cube
    min_neigh_dist = sqrt(3)*median(min_neigh_dist);
    % Iterations
    adja_restr = true;
    new_assignment = old_assignment;
    for iter = 1:max_iter
        
    %                 Cluster Processing Order Definition                %
    
        proc_order = zeros(1,num_clust);
        clus_occupation = hist(new_assignment,num_clust);
        % In the second iteration, determine the smallest cluster(s) for
        % MIC or SIC cases
        if extra > 0 && iter == 2
            [~, idx_smaller_cluster] = sort(clus_occupation);
            idx_smaller_cluster(num_smaller + 1:end) = [];
        end
        % Check if iteration number is even or odd, to select the first cluster as min or max
        if mod(iter, 2) == 0
            % If smaller cluster(s), in MIC or SIC cases, is (are) already full,
            % ignore it (them)
            if sum(clus_occupation(idx_smaller_cluster) == min_size_cluster)
                clus_occupation(idx_smaller_cluster(clus_occupation(idx_smaller_cluster) == min_size_cluster)) = max(clus_occupation);
            end
            % Choose the cluster with least points in current iteration as
            % the first to process
            [~, proc_order(1)] = min(clus_occupation);
        else
            % If the most populated cluster are already at the maximum size
            % (full), ignore them
            if max(clus_occupation) == max_size_cluster
                clus_occupation(clus_occupation == max_size_cluster) = 0;
            end
            % Choose the cluster with most points in current iteration as
            % the first to process
            [~, proc_order(1)] = max(clus_occupation);
        end
        % Sort remaining clusters by distance to the first cluster
        [~,current_dist_idx] = pdist2(centroids,centroids(proc_order(1),:),'euclidean','Smallest',num_clust);
        aux = setdiff(current_dist_idx, proc_order(1), 'stable');
        % The processing order corresponds to the sorted clusters
        proc_order(2:num_clust) = aux;

    %                 Cluster-to-Cluster Point Transfer                  %
    
        % Clear full clusters list
        full_clusters = [];
        % Process each cluster in order
        for i = 1:num_clust
            % Check the target size for the current cluster
            if ismember(proc_order(i), idx_smaller_cluster)
                size_cluster = min_size_cluster;
            else
                size_cluster = max_size_cluster;
            end
            % Get coordinates for current cluster
            cluster = coordinates(new_assignment == proc_order(i), :);
            % Keep list of IDs for points in current cluster
            assign_id = find(new_assignment == proc_order(i));
            % Check if cluster has more points than the target
            if length(cluster) > size_cluster                
                excessive_points = length(cluster) - size_cluster;
                % Determine neighbors eligible for transfers %
                % Get the nearest neighbors clusters
                [~,neighb_idx] = pdist2(centroids,centroids(proc_order(i),:),'euclidean','Smallest',num_nn + 1); neighb_idx(1) = [];
                % Check if adjacency restriction is still in place
                if adja_restr
                    % Add current cluster to adjacent neighborhood
                    temp_neighborhood = cluster;
                    neighbor_list = [];
                    % For each nearest neighbor, check if it is adjacent to
                    % the current cluster, or adjacent to another adjacent
                    % cluster (i.e., adjacent to the neighborhood)
                    for j=1:length(neighb_idx)
                        % Get points for current nearest neighbor
                        next_neighb = coordinates(new_assignment == neighb_idx(j), :);
                        % Determine closest distances between all points of
                        % nearest neighbor and the current neighborhood
                        [dist_neigh,~] = pdist2(temp_neighborhood,next_neighb,'euclidean','Smallest',1);
                        % If the minimum distance is below or equal the
                        % threshold, it is considered adjacent
                        if min(dist_neigh) <= min_neigh_dist
                            % Add cluster to the adjacent neighborhood
                            temp_neighborhood = [temp_neighborhood; next_neighb];
                            neighbor_list = [neighbor_list; neighb_idx(j)];
                        end
                    end
                else
                    % If the adjacency restriction was lifted, consider all
                    % nearest neighbors
                    neighbor_list = neighb_idx;
                end
                % Discard full clusters already processed in this iteration
                neighbor_list = setdiff(neighbor_list, full_clusters, 'stable');
                % If there are no eligible clusters, continue to the next
                % cluster
                if isempty(neighbor_list)
                    continue;
                end
                % Transfer points to neighboring clusters %
                % Compute distances from every point of the cluster to the eligible clusters centroids
                [p_out_dist,p_out_dist_idx] = pdist2(centroids(neighbor_list,:),cluster,'euclidean','Smallest',length(neighbor_list));
                % Flatten matrices
                assign_id = reshape(repmat(assign_id',length(neighbor_list),1),1,[]); p_out_dist = reshape(p_out_dist,1,[]); p_out_dist_idx = reshape(p_out_dist_idx,1,[]);
                % Sort all possible transfers by shortest distance
                [~,points_to_move] = sort(p_out_dist,'ascend');
                k = 0;
                already_moved = [];
                changed = [];
                % Transfer points out of the current cluster until the
                % cluster is full
                for j = 1:excessive_points
                    % Run through possible transfers
                    while(k < length(points_to_move))
                        k = k + 1;
                        % Get cluster ID of current possible transfer
                        new_id = neighbor_list(p_out_dist_idx(points_to_move(k)));
                        % Check if the transfer is possible (neighbor
                        % cluster not full, and current point to move has
                        % not already moved to another cluster)
                        if ~ismember(new_id, full_clusters) && ~ismember(assign_id(points_to_move(k)), already_moved)
                            % Perform point transfer
                            new_assignment(assign_id(points_to_move(k))) = new_id;
                            % Save list of points that have been moved
                            already_moved = [already_moved assign_id(points_to_move(k))];
                            % Save list of clusters that changed
                            changed = [changed new_id];
                            break;
                        end
                    end
                end
                % Update centroid coordinates for changed clusters
                changed = unique(changed);
                for j = 1:length(changed)
                    centroids(changed(j),:) = mean(coordinates(new_assignment == changed(j), :));
                end
            % Check if cluster has fewer points than the target
            elseif length(cluster) < size_cluster
                lacking_points = abs(length(cluster) - size_cluster);
                % Determine neighbors eligible for transfers %
                % Get the nearest neighbors clusters
                [~,neighb_idx] = pdist2(centroids,centroids(proc_order(i),:),'euclidean','Smallest',num_nn + 1); neighb_idx(1) = [];
                % Check if adjacency restriction is still in place
                if adja_restr
                    % Add current cluster to adjacent neighborhood
                    temp_neighborhood = cluster;
                    neighbor_list = [];
                    % For each nearest neighbor, check if it is adjacent to
                    % the current cluster, or adjacent to another adjacent
                    % cluster (i.e., adjacent to the neighborhood)
                    for j = 1:length(neighb_idx)
                        % Get points for current nearest neighbor
                        next_neighb = coordinates(new_assignment == neighb_idx(j), :);
                        % Determine closest distances between all points of
                        % nearest neighbor and the current neighborhood
                        [dist_neigh,~] = pdist2(temp_neighborhood,next_neighb,'euclidean','Smallest',1);
                        % If the minimum distance is below or equal the
                        % threshold, it is considered adjacent
                        if min(dist_neigh) <= min_neigh_dist
                            % Add cluster to the adjacent neighborhood
                            temp_neighborhood = [temp_neighborhood; next_neighb];
                            neighbor_list = [neighbor_list; neighb_idx(j)];
                        end
                    end
                else
                    % If the adjacency restriction was lifted, consider all
                    % nearest neighbors
                    neighbor_list = neighb_idx;
                end
                % Discard full clusters already processed in this iteration
                avail_neighb = setdiff(neighbor_list, full_clusters, 'stable');
                % If there are no eligible clusters, continue to the next
                % cluster
                if isempty(avail_neighb)
                    continue;
                end
                % Get points from neighbouring clusters %
                % Get all points from the eligible clusters
                cand_idx = find(sum(new_assignment == avail_neighb',2) == 1);
                candidates = coordinates(cand_idx, :);
                % Order candidates by distance to the current cluster centroid
                [~,p_in_dist_idx] = pdist2(candidates,centroids(proc_order(i),:),'euclidean','Smallest',lacking_points);
                changed = proc_order(i);
                % Transfer points into the current cluster until the
                % cluster is full
                for j = 1:length(p_in_dist_idx)
                    % Only transfer if the neighboring cluster where the possible
                    % transfer point is from has at least 10% of the target
                    % PPC, to avoid emptying it
                    if sum(new_assignment == new_assignment(cand_idx(p_in_dist_idx(j)))) > round(0.1*ppc)
                        % Save list of clusters that changed
                        changed = [changed new_assignment(cand_idx(p_in_dist_idx(j)))];
                        % Perform point transfer
                        new_assignment(cand_idx(p_in_dist_idx(j))) = proc_order(i);
                    end
                end
                % Update centroid coordinates for changed clusters
                changed = unique(changed);
                for j=1:length(changed)
                    centroids(changed(j),:) = mean(coordinates(new_assignment == changed(j), :));
                end
            end
            % Check if the current cluster achieved its target size (with
            % or without transfers); If so, flag it as a full cluster
            if length(find(new_assignment == proc_order(i))) == size_cluster
                full_clusters = [full_clusters proc_order(i)];
            end
        end
        % Track progress of algorithm for debugging purposes
        if debug
            for j=1:num_clust
                cluster = coordinates(new_assignment == j, :);
                clust_disp(j,:,iter+1) = var(cluster);
                clust_bb(j,:,iter+1) = range(cluster);
            end
        end

    %                       Check STOP condition                         %
    
        progress(:,iter+1) = hist(new_assignment,num_clust)';
        if (max(progress(:,iter+1)) == max_size_cluster && min(progress(:,iter+1)) == min_size_cluster)
            break;
        end        

    %          Check if adjacency condition needs to be removed          %
    
        % Check if the total number of "misplaced" points did not decrease
        % from the previous iteration. "Misplaced" points are the points
        % above or below the target cluster size.
        if (sum(abs(progress(:,iter+1)-ppc)) >= sum(abs(progress(:,iter)-ppc)))
            if (adja_restr == true)
                temp_calc = sum(abs(progress - ppc));
                temp_idx = find(temp_calc == temp_calc(iter+1), 1, 'first');
                % Check if the current configuration of number of points in
                % each cluster has been repeated; or if the number of
                % "misplaced" points has not improved in over 5 iterations.
                % This indicates if the algorithm has stalled
                if ~isequal(unique(progress(:,temp_idx:iter+1)','rows','stable'), progress(:,temp_idx:iter+1)') || (iter+1 - temp_idx > 5)
                    if debug
                        % Record the iteration number in which it was
                        % deactivated
                        iter_adj_rest = iter;
                    end
                    % Remove adjacency condition
                    adja_restr = false;
                end
            else
                % If the adjacency condition is not active anymore,
                % increase the number of nearest neighbors to check
                num_nn = num_nn + 10;
            end
        end
    end
    progress(:,iter+2:end) = [];
    if debug
        fprintf("Number of iterations in iterative phase: %d\n",iter);
        clust_disp(:,:,iter+2+2:end) = [];
        clust_bb(:,:,iter+2+2:end) = [];
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        Clustering Refinement                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    ref_assignment = new_assignment;
    if debug
        % Measure cluster dispersion before refinement
        sum_dists = zeros(101,1);
        temp_sum = zeros(num_clust,1);
        for j = 1:num_clust
            cluster = coordinates(ref_assignment == j, :);
            point_dist = pdist2(mean(cluster),cluster,'euclidean');
            temp_sum(j) = sum(point_dist);
        end
        sum_dists(1) = mean(temp_sum);
    end
    % Continue iterating for all clusters until no swaps are performed
    for iter = 1:100
        change_flag = false;
        % Process each cluster
        for i = 1:num_clust
            % Get the 10 nearest cluster centroids
            [~,aaa_idx] = pdist2(centroids,centroids(i,:),'euclidean','Smallest',11); aaa_idx(1) = [];
            % For each nearest neighbor
            for j = 1:length(aaa_idx)
                % Get points of current cluster A
                assign_id_A = find(ref_assignment == i);
                clusterA = coordinates(ref_assignment == i, :);
                % Get points of neighboring cluster B
                assign_id_B = find(ref_assignment == aaa_idx(j));
                clusterB = coordinates(ref_assignment == aaa_idx(j), :);
                % Compute distance of points in A to centroid of A
                cAA_dist = pdist2(clusterA,centroids(i,:),'euclidean');
                % Compute distance of points in A to centroid of B
                cAB_dist = pdist2(clusterA,centroids(aaa_idx(j),:),'euclidean');
                % Compute distance of points in B to centroid of A
                cBA_dist = pdist2(clusterB,centroids(i,:),'euclidean');
                % Compute distance of points in B to centroid of B
                cBB_dist = pdist2(clusterB,centroids(aaa_idx(j),:),'euclidean');
                % Get and sort candidate points of A by preference to be in B
                [auxAB,candidatesAB] = sort(cAA_dist - cAB_dist, 'descend');
                % Get and sort candidate points of B by preference to be in A
                [auxBA,candidatesBA] = sort(cBB_dist - cBA_dist, 'descend');
                % Run through both candidates lists until one of them is
                % over, or until neither candidate prefers to change
                for k = 1:min(max(find(auxAB <= 0, 1) - 1,find(auxBA <= 0, 1) - 1), min(length(auxAB), length(auxBA)))
                    % Do not swap if: (a) only one point wants to change, but the benefit of
                    % changing it is lower than the loss of changing the
                    % other point; (b) neither point wants to change
                    if ((auxAB(k) < 0 && abs(auxAB(k)) >= auxBA(k)) ||...
                            (auxBA(k) < 0 && abs(auxBA(k)) >= auxAB(k)) ||...
                            (auxAB(k) < 0 && auxBA(k) < 0))
                        break
                    end
                    % Swap points
                    ref_assignment(assign_id_A(candidatesAB(k))) = aaa_idx(j);
                    ref_assignment(assign_id_B(candidatesBA(k))) = i;
                    change_flag = true;
                end
                % Update centroid coordinates
                centroids(i,:) = mean(coordinates(ref_assignment == i, :));
                centroids(aaa_idx(j),:) = mean(coordinates(ref_assignment == aaa_idx(j), :));
            end                
        end
        if debug
            % Measure cluster dispersion in each iteration
            for j=1:num_clust
                cluster = coordinates(ref_assignment == j, :);
                point_dist = pdist2(mean(cluster),cluster,'euclidean');
                temp_sum(j)=sum(point_dist);
            end
            sum_dists(iter+1) = mean(temp_sum);
        end
        % Terminate refinement if no swaps occurred
        if ~change_flag
            break
        end
    end
    if debug
        fprintf("Number of iterations in refinement phase: %d\n",iter);
        for j=1:num_clust
            cluster = coordinates(ref_assignment == j, :);
            clust_disp(j,:,end-1) = var(cluster);
            clust_bb(j,:,end-1) = range(cluster);
        end
        sum_dists(iter+2:end) = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Cluster Padding                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    if padding
        % Get clusters with fewer points than the target size
        idx = find(progress(:,end) < ppc);
        % For each of the smaller clusters
        for i=1:length(idx)
            % Compute the number of points missing
            extra = ppc - progress(idx(i),end);
            if(extra < ppc)
                % Get cluster points
                cluster = coordinates(ref_assignment == idx(i), :);
                % Randomly repeat points until cluster is full
                while(extra > 0)
                    rem_idx = randsample(1:progress(idx(i),end), min(extra, progress(idx(i),end)));
                    coordinates = [coordinates; cluster(rem_idx, :)];
                    ref_assignment = [ref_assignment; idx(i).*ones(length(rem_idx),1)];
                    extra = extra - length(rem_idx);
                end
            end
        end
        if debug
            for j=1:num_clust
                cluster = coordinates(ref_assignment == j, :);
                clust_disp(j,:,end) = var(cluster);
                clust_bb(j,:,end) = range(cluster);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Write Results                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % Create folder for final clusters
    if(~isfolder(outfolder + 'Refined'))
        mkdir(outfolder + 'Refined');
    end
    if debug
        % Create folder also for clusters after initial and iterative phases
        if(~isfolder(outfolder + 'Iterative'))
            mkdir(outfolder + 'Iterative');
        end
        if(~isfolder(outfolder + 'Initial'))
            mkdir(outfolder + 'Initial');
        end
        % Output additional results
        if nargout >= 1
            varargout{1} = progress;
        end
        if nargout >= 2
            varargout{2} = clust_disp;
        end
        if nargout >= 3
            varargout{3} = clust_bb;
        end
        if nargout >= 4
            varargout{4} = sum_dists;
        end
        if nargout >= 5
            varargout{5} = iter_adj_rest;
        end
        if nargout >= 6
            for ii=6:nargout
                varargout{ii} = [];
            end
        end
    else
        if nargout >= 1
            varargout{1} = progress;
        end
        if nargout >= 2
            for ii=2:nargout
                varargout{ii} = [];
            end
        end
    end
    for j=1:num_clust
        % Randomly choose a color for each cluster
        col_rgb = randi(256, 1, 3) - 1;
        % Get cluster points, and write it to a PLY file
        cluster = coordinates(ref_assignment == proc_order(j), :);
        pcwrite(pointCloud(cluster, 'Color', repmat(uint8(col_rgb),length(cluster),1)), sprintf('%s/%s_%05d.ply', outfolder + 'Refined',outpcname,j),'PLYFormat','binary');
        % If debug is activated, also write intermediate clusters
        if debug
            cluster = coordinates(new_assignment == proc_order(j), :);
            pcwrite(pointCloud(cluster, 'Color', repmat(uint8(col_rgb),length(cluster),1)), sprintf('%s/%s_%05d.ply', outfolder + 'Iterative',outpcname,j),'PLYFormat','binary');
            cluster = coordinates(old_assignment == proc_order(j), :);
            pcwrite(pointCloud(cluster, 'Color', repmat(uint8(col_rgb),size(cluster,1),1)), sprintf('%s/%s_%05d.ply', outfolder + 'Initial',outpcname,j),'PLYFormat','binary');
        end
    end
end
