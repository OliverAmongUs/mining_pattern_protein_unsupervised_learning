 function [ match_score_matrix ] = graph_matching_new( ARG1, ARG2, train )
 %  We use the structure motif learning for rigid objects algorithms
 %  In this function we calculate the sum of pair-wise matching score
 %  for each sample-component matching
 %  Also, we meed to record 
 %  1. the matchings between the chosen protein graph 
 %  and all other protein graphs
 %  2. The transformation matrixes between the chosen protein graph and all
 %  other protein graphs
    
        
    
    % number of nodes in each graph
    size_1 = ARG1.num_nodes;
    size_2 = ARG2.num_nodes;
    
    % initialize two cells to store subgraphs, 
    % subgraph_1 for ARG1 and subgraph_2 for ARG2
    subgraph_1 = cell(size_1,1);
    subgraph_2 = cell(size_2,1);
    
    % How to initialize properly(diagonal should be zero)
    % Check if this is correct
    match_score_matrix = ones([size_1, size_2])/(size_1*size_2);
    
    % Matching
    for i = 1:ARG1.num_nodes
        for j = 1:ARG2.num_nodes
            
            % Find the subgraph centered at 
            % ith node(in ARG1) and jth node(in ARG2)
            
            match_score_local = graph_matching_subgraph(subgraph_1, subgraph_2, ARG1, ARG2, i, j);
            
            % make sure this is correct
            match_score_matrix(i,j) = match_score_local
            
        end
    end
    
     function [match_score_local] = graph_matching_subgraph(subgraph_1, subgraph_2,ARG1, ARG2, i,j)
         % First we focus on ith (jth) node and its edge in ARG1(ARG2)
         subgraph_1_tmp = ARG1.edges_matrix(:,i);
         subgraph_2_tmp = ARG2.edges_matrix(:,j);
         
         % Then we extract the subgraph: we first find all the nodes that
         % are connected to i(j), and then extract these nodes together
         % with their edges (edges that both nodes are in this subgrah)
         
         % We set the subgraph_1_tmp[i,i]([j,j])to be one because we want to
         % include itself into the subgraph
         
         subgraph_1_tmp(i) = 1;
         subgraph_2_tmp(j) = 1;
         
         % https://www.mathworks.com/matlabcentral/answers/179176-how-can-i-find-the-rows-with-l-non-zero-elements-in-a-matrix-vital
         rowIdcs1 = (find(subgraph_1_tmp~=0)).';
         rowIdcs2 = (find(subgraph_2_tmp~=0)).';
        
         % Two subgraphs extracted centered at i and j
         subgraph_1{i} = ARG1.edges_matrix(rowIdcs1,rowIdcs1);
         subgraph_2{j} = ARG2.edges_matrix(rowIdcs2,rowIdcs2);
         
         subgraph_nodes_1 = num2cell(ARG1.nodes_vector(rowIdcs1).');
         
         %!!!!This sentence is not correct
         subgraph_nodes_2 = num2cell(ARG2.nodes_aa_index(rowIdcs2).');
         
         a = ARG(num2cell(subgraph_1{i}),subgraph_nodes_1);
         b = mdl_ARG(ARG(num2cell(subgraph_2{j}),subgraph_nodes_2));
         %match_score_local = prctile(graph_matching(a,b,train),50)
         match_score_local = median(reshape(match_score(a,b),1,[]))
         
         
     end
     function [matching_score] = match_score(ARG1,ARG2)
         % show procedure flag
         sFlag = 0;
         % set up condition and variable
         % beta is the converging for getting the maximize number
         beta_0 = 0.1;
         beta_f = 20;   % original is 10
         beta_r = 1.025; % original is 1.075

         % I control the iteration number for each round
         I_0 = 20;  % original is 4
         I_1 = 200;   % original is 30

         % e control a range
         e_B = 0.5;  % original is 0.5
         e_C = 0.05;   % original is 0.05

         % e_cov to handle singularity
         e_cov = 0.01;

         % node attriubute compatability weight
         alpha = 0.1; % edge
         delta = 1; % node
         gamma = 1; %linear

         % the size of the real matchin matrix
         A=ARG1.num_nodes;
         I=ARG2.num_nodes;
         real_size = [A,I];
         
         % adjust ARG2
         M = ARG2.edges_matrix(1:I,1:I,:);
    
         % covariance of ARG2
         M_C = ARG2.edges_cov(1:I,1:I,:);
    
         % eliminate the extra information
         V = ARG2.nodes_vector(1:I,:);
         % the size of the matrix with slacks
         augment_size = real_size+[1,1];

         % initial beta to beta_0
         beta = beta_0;

         % nil node compatibility percentage
         prct = 105;

         % stochastic level
         s_level = 1;

         % pre-calculate the node compatability
         C_n=zeros(A+1,I+1);

         % create an function handle for calculating compatibility
         % calculate the compatibility
         for a = 1:A
             for i = 1:I
                 C_n(a,i) = node_compatibility(ARG1.nodes_vector(a,:),V(i,:));
             end
         end

         % calculate nil compatibility
         C_n(A+1, 1:I)=max(0,prctile(C_n(1:A,1:I),min(prct,100),1)*max(prct/100,1));
         C_n(1:A, I+1)=max(0,prctile(C_n(1:A,1:I),min(prct,100),2)*max(prct/100,1));
         C_n(A+1, I+1)=0;
         
         % How do we define a matching score
         if sFlag
             figure()
         end
         
         %start matching
         while beta<beta_f %do A until beta is less than beta_f
             
             converge_B = 0; % a flag for terminating process B
             I_B = 0;        % counting the iteration of B
             
             while ~converge_B && I_B <= I_0 % do B until B is converge or iteration exceeds
                 
                 
    function [score] = node_compatibility(atr1, atr2)
        score = atr2(atr1);
    end
 end