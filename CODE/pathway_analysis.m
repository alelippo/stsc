% MATLAB script used to generate the data in Section 6 of the paper:
%
%       "Source-to-sink communicability: a new centrality measure for
%        directed acyclic networks", submitted.
%
%  Requirements: "DATASET" folder, "stsc.m" function  
%   
%  Daniele Bertaccini, Luigi Chiricosta and Alessandro Filippo, 2025.

clear; clc; close all; lastwarn(''); 

% Choice of the pathway

name = 'hsa05022';
%name = 'hsa05200';
%name = 'hsa05010';
%name = 'hsa04010';

disp(name)
filename = ['adjacency_matrix_',name,'_R.tsv'];
data = readtable(filename, "FileType","text",'Delimiter', '\t', 'VariableNamingRule','preserve');
n = size(data,2);

% Extract node names
%NOTE: we have to react if node names are longer than 63 char (MATLAB LIMIT)  
[~, warnId] = lastwarn;
if isempty(warnId) % if there are no warnings
    nodeNames = data.Properties.VariableNames;
else
    nodeNames = data.Properties.VariableDescriptions;
end
nodeNames = nodeNames(2:n);

% Table to matrix conversion
A = sparse(data{:,2:n});

% Create the graph object
G = digraph(A, nodeNames);

% Remove isolated nodes 
nodeDegrees = indegree(G) + outdegree(G);
isolatedNodes = find(nodeDegrees == 0);
G = rmnode(G, isolatedNodes);
nodeNames = G.Nodes.Name;
A = adjacency(G);
n = size(A,1);

% How to deal with cyclic digraphs
% NOTE: The graph is now topologically sorted, i.e., A is upper triangular.
if hascycles(G)
    disp('--')
    disp('The graph is not ACYCLIC!')
    disp('Moving to the condensated digraph')
    disp('--')
    C = condensation(G); 
    bins = conncomp(G);
    temp = accumarray(bins.', (1:length(bins)).', [], @(IDX) {G.Nodes.Name(IDX)});
    nodeNames = cellfun(@(CS) strjoin(CS, ';'), temp, 'uniform', 0);
    A = adjacency(C);
    n = size(A,1);
    C.Nodes.Name = nodeNames;
    G = C;
end

% Graph info
m = G.numedges;
avg_d = m/n;
diam = my_diameter(G);

% Number of "top nodes" to display 
k = 5;

% DEGREE

ind = sum(A,1)'; % column sum
outd = sum(A,2);  % row sum
[values,I_deg] = maxk(ind+outd,n);
values_deg = values(1:k);
top_deg_sum = nodeNames(I_deg(1:k));

% KATZ 

alpha = 0.85;
I = eye(n);
b = ones(n,1);
S = I-alpha*A;
inx = S'\b;
outx = S\b;
[values,I_katz] = maxk(inx + outx,n);
values_katz = values(1:k);
top_katz_sum = nodeNames(I_katz(1:k));

% PAGERANK and REVERSED PAGERANK

Gt = digraph(A');
pi = centrality(G,'pagerank');
pi_rev = centrality(Gt,'pagerank');
[values,I_pgr] = maxk(pi + pi_rev,n);
values_pgr = values(1:k);
top_pgr_sum = nodeNames(I_pgr(1:k));

% HUB AND AUTHORITIES

h = centrality(G,'hubs');
a = centrality(G,'authorities');
[values,I_hub_auth] = maxk(h + a,n);
values_hub_auth = values(1:k);
top_hub_auth_sum = nodeNames(I_hub_auth(1:k));

% STS communicability

x = stsc(A);
[values,I_sts] = maxk(x,n);
values_sts = values(1:k);
top_sts = nodeNames(I_sts(1:k));

disp(table(n,m,avg_d,diam))
disp(table(top_sts,values_sts)) 

disp('--')

% KENDAL COEFFICIENT

disp('Kendall coefficients')
disp([' deg = ',num2str(corr(I_deg,I_sts,'type','Kendall'))])
disp([' pgr = ',num2str(corr(I_pgr,I_sts,'type','Kendall'))])
disp([' katz = ',num2str(corr(I_katz,I_sts,'type','Kendall'))])
disp([' hub_auth = ',num2str(corr(I_hub_auth,I_sts,'type','Kendall'))])
disp('--')

% INTERSECTION SIMILARITY

disp('Top k intersection similarity')
disp([' deg = ',num2str(my_top_k_intersection_similarity(I_sts(1:k),I_deg(1:k)))])
disp([' pgr = ',num2str(my_top_k_intersection_similarity(I_sts(1:k),I_pgr(1:k)))])
disp([' katz = ',num2str(my_top_k_intersection_similarity(I_sts(1:k),I_katz(1:k)))])
disp([' hub_auth = ',num2str(my_top_k_intersection_similarity(I_sts(1:k),I_hub_auth(1:k)))])
disp('--')

function I = my_top_k_intersection_similarity(x,y,k)
% Function that computes the "top k Intersection similarity" of two ordered
% lists as described in R. Fagin, R. Kumar, D. Sivakumar. Comparing top k 
% lists. SIAM J. Discret. Math., 17:134–160, 2003
    if nargin<3
        k = min(length(x),length(y));
    end
    I = 0;     
    for j = 1:k
        I = I + length(setxor(x(1:j),y(1:j)))/(2*j); 
    end    
    I = I/k; 
end

