# yaoxiaojing

Extracting useful spatial co-location patterns from urban service facilities can help planners allocate limited resources effectively. These facilities are mostly distributed within man-made spatial fields with road-network constraints. To promote urban-space adaptivity, some current co-location algorithms for this network circumstance have designed with distance decay effects and topological relationships of roads. However, they neglect the traffic direction, which affects the result accuracy. Moreover, the efficiency problem is more severe compared with the traditional algorithm (i.e., no constraints). To address these problems, we propose an efficient maximal co-location mining algorithm, with directed road-network constraints and spatial-continuity consideration (CMDS). To improve the accuracy, we design a network-based prevalence index, combined with both distance decay effects and road direction interference, to measure the significance of a pattern. To promote the execution speed, we use a key-node-separating approach and an improved shortest-path batch task for the co-location mining process. The experiments with both the synthetic and real datasets show that the CMDS algorithm is more efficient and accurate than the state-of-the-art network co-location when applied to problems in an urban space. 


/* the useage of the files

1. The core algorithm is coded by Matlab with a file called "net_kde_colocation_4.m".

2. "Tree/@tree" is a package to deal with tree-based implementations.

3. "网络test" contains experimental datasets.
