function x = stsc(A)
%%STSC  computes the vector x such that x_w is the source-to-sink 
%       communicability of node w, i.e., 
%       x = s^T (I-A)^{-1} ∘ (I-A)^{-1} t
%       requiring A being the adjacency matrix of an acyclic digraph G.
%
%   Reference: ####
%   
%   Daniele Bertaccini, Luigi Chiricosta and Alessandro Filippo, 2025.

s = sum(A,1)' == 0;             % source nodes
t = sum(A,2) == 0;              % sink nodes
M = speye(size(A,1)) - A;   

if hascycles(digraph(A))        % Check if the graph is acyclic
    error('The digraph is not ACYCLIC.')
end

y = M'\s;                       % Solving the linear systems
z = M\t;

x = (y.*z) / ((s')*z) ;         % stsc computation

end