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
    subgraph_1 = cell(size_1);
    subgraph_2 = cell(size_2);
    
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
        
         % Two subgraphs extrated centered at i and j
         subgraph_1{i} = ARG1.edges_matrix(rowIdcs1,rowIdcs1);
         subgraph_2{j} = ARG2.edges_matrix(rowIdcs2,rowIdcs2);
         
         subgraph_nodes_1 = num2cell(ARG1.nodes_vector(rowIdcs1).')
         subgraph_nodes_2 = num2cell(ARG2.nodes_vector(rowIdcs2).')
         
         a = ARG(num2cell(subgraph_1{i}),subgraph_nodes_1);
         b = mdl_ARG(ARG(num2cell(subgraph_2{j}),subgraph_nodes_2));
         %match_score_local = prctile(graph_matching(a,b,train),50)
         match_score_local = median(reshape(graph_matching(a,b,train),1,[]))
         
         
     end
 end
 
         

         
         
         
         
         
         
         
         
            
            
     
            
       
            
            
            
            
            
            
            
            
    