 function [ match_matrix ] = graph_matching_new( ARG1, ARG2, train )
    
    % number of nodes in each graph
    size_1 = ARG1.num_nodes;
    size_2 = ARG2.num_nodes;
    
    % initialize two cells to store subgraphs, 
    % subgraph_1 for ARG1 and subgraph_2 for ARG2
    subgraph_1 = cell(size_1);
    subgraph_2 = cell(size_2);
    
    % How to initialize properly(diagonal should be zero)
    match_score_matrix = ones([size_1, size_2])/(size_1*size_2);
    
    % Matching
    for i = 1:ARG1.num_nodes
        for j = 1:ARG2.num_nodes
            
            % Find the subgraph centered at 
            % ith node(in ARG1) and jth node(in ARG2)
            
            match_score_local = graph_matching_subgraph(ARG1, ARG2, i, j);
            
            %null-node-score?
        end
    end
    
     function [match_score_local] = graph_matching_subgraph(ARG1, ARG2, i,j)
         % First we focus on ith (jth) node and its edge in ARG1(ARG2)
         subgraph_1_tmp = ARG1.edges_matrix(:,i);
         subgraph_2_tmp = ARG2.edges_matrix(:,j);
         
         % Then we extract the subgraph: we first find all the nodes that
         % are connected to i(j), and then extract these nodes together
         % with their edges (edges that both nodes are in this subgrah)
         
         % We set the subgraph_1_tmp[i,i]([j,j])to be one because we want to
         % include itself into the subgraph
         
         subgraph_1_tmp(i,i) = 1;
         subgraph_2_tmp(j,j) = 1;
         
         % https://www.mathworks.com/matlabcentral/answers/179176-how-can-i-find-the-rows-with-l-non-zero-elements-in-a-matrix-vital
         rowIdcs1 = (find(subgraph_1_tmp~=0)).';
         rowIdcs2 = (find(subgraph_2_tmp~=0)).';
         
         
         subgraph_1{i} = ARG1.edges_matrix(rowIdcs1,rowIdcs1);
         subgraph_2{i} = ARG1.edges_matrix(rowIdcs2,rowIdcs2);
     end
 end
 
         

         
         
         
         
         
         
         
         
            
            
     
            
       
            
            
            
            
            
            
            
            
    