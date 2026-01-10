[INFORMATION] current ppt index's format: Table
[LOAD GRAPH] start!!! graph: example filepath: ../datasets/example.txt
graph size: 12
edge size: 18
[LOAD GRAPH] success!!!
[CONSTRUCT INDEX] start building BOTTOM-K-SKETCH
[CONSTRUCT INDEX] finish building BOTTOM-K-SKETCH
[TIME COST] 2.8007e-05 s.
[CONSTRUCT INDEX] start building Neighbor Order
[CONSTRUCT INDEX] finish building Neighbor Order
[TIME COST] 1.8006e-05 s.
[CONSTRUCT Index] start PPT decomposition
[CONSTRUCT INDEX] finish PPT decomposition
[TIME COST] 3.8111e-05 s.
[CONSTRUCT INDEX] start building componets
[CONSTRUCT INDEX] finish building componets
[TIME COST] 8.6849e-05 s.
[PRINT NODE NUM] total point (leaf node): 12
[PRINT TIME] total building time: 0.000170973
[PRINT SPACE INFORMATION] core vertices space cost: 9.15527e-05.MB
[PRINT SPACE INFORMATION] non core vertices space cost: 4.57764e-05.MB
[PRINT SPACE INFORMATION] edges space cost: 9.15527e-05.MB
[PRINT SPACE INFORMATION] connectivity space cost: 4.57764e-05.MB
[PRINT PPT SPACE]: total space cost: 0.000411987.MB
[PRINT] print NEIGHBOR ORDER start !!!
0: (0.666667,2)(0.6,3)(0.6,1)
1: (0.666667,2)(0.6,0)(0.5,4)
2: (0.666667,0)(0.666667,1)(0.666667,3)(0.571429,5)(0.571429,4)
3: (0.666667,2)(0.6,0)(0.5,5)
4: (0.571429,2)(0.5,1)(0.428571,5)(0.4,6)
5: (0.571429,2)(0.5,3)(0.428571,4)(0.333333,7)
6: (0.4,4)
7: (0.4,8)(0.333333,5)
8: (0.6,10)(0.6,9)(0.4,7)
9: (1,10)(0.75,11)(0.6,8)
10: (1,9)(0.75,11)(0.6,8)
11: (0.75,9)(0.75,10)
[PRINT] print NEIGHBOR ORDER finish !!!
[PRINT] now print the similarity interval
[0,0.1) [0.1,0.2) [0.2,0.3) [0.3,0.4) [0.4,0.5) [0.5,0.6) [0.6,0.7) [0.7,0.8) [0.8,0.9) [0.9,1) 
[PRINT] print finish !!!
[PRINT] print INTERVAL INDEX start !!!
0: 0 0 0 0 0 0 3 0 0 0 
1: 0 0 0 0 0 1 2 0 0 0 
2: 0 0 0 0 0 2 3 0 0 0 
3: 0 0 0 0 0 1 2 0 0 0 
4: 0 0 0 0 2 2 0 0 0 0 
5: 0 0 0 1 1 2 0 0 0 0 
6: 0 0 0 0 1 0 0 0 0 0 
7: 0 0 0 1 1 0 0 0 0 0 
8: 0 0 0 0 1 0 2 0 0 0 
9: 0 0 0 0 0 0 1 1 0 1 
10: 0 0 0 0 0 0 1 1 0 1 
11: 0 0 0 0 0 0 0 2 0 0 
[PRINT] print INTERVAL INDEX finish !!!
[PRINT] now print the map between PPT and vertex !!!
All PPT: 
(6,4) vertices: 0 2 9 10 
(6,3) vertices: 1 3 8 
(5,4) vertices: 1 3 
(5,6) vertices: 2 
(5,3) vertices: 4 5 
(4,5) vertices: 4 
(4,4) vertices: 5 8 
(3,5) vertices: 5 
(4,2) vertices: 6 7 
(3,3) vertices: 7 
(9,2) vertices: 9 10 
(7,3) vertices: 9 10 11 
max mu under each epsilon: 
epsilon bucket 0 : 1
epsilon bucket 1 : 1
epsilon bucket 2 : 1
epsilon bucket 3 : 5
epsilon bucket 4 : 5
epsilon bucket 5 : 6
epsilon bucket 6 : 4
epsilon bucket 7 : 3
epsilon bucket 8 : 1
epsilon bucket 9 : 2
PPT matrix: 



-1 9 -1 7 
8 -1 6 5 
-1 4 2 -1 3 
-1 1 0 
-1 11 

10 
Vertex 2 PPT: 
vertex: 0 PPT: (6,4)
vertex: 1 PPT: (6,3)(5,4)
vertex: 2 PPT: (6,4)(5,6)
vertex: 3 PPT: (6,3)(5,4)
vertex: 4 PPT: (5,3)(4,5)
vertex: 5 PPT: (5,3)(4,4)(3,5)
vertex: 6 PPT: (4,2)
vertex: 7 PPT: (4,2)(3,3)
vertex: 8 PPT: (6,3)(4,4)
vertex: 9 PPT: (9,2)(7,3)(6,4)
vertex: 10 PPT: (9,2)(7,3)(6,4)
vertex: 11 PPT: (7,3)
[PRINT] print finish !!!
[RPINT] now print all ppt after construct componets
=================================================================================
id: 0 epsilon: 6 mu: 4
vertices: 0 2 9 10 
subgraph: (cores: 9 10 non-cores: 11 8 )(cores: 0 2 non-cores: 3 1 )
connectivity: 
=================================================================================
id: 1 epsilon: 6 mu: 3
vertices: 1 3 8 
subgraph: (cores: 8 non-cores: )(cores: 1 3 non-cores: )
connectivity: (1,2):0.666667 (8,10):0.6 
=================================================================================
id: 2 epsilon: 5 mu: 4
vertices: 1 3 
subgraph: (cores: 1 3 non-cores: )
connectivity: (1,2):0.666667 
=================================================================================
id: 3 epsilon: 5 mu: 6
vertices: 2 
subgraph: (cores: 2 non-cores: 0 1 3 5 4 )
connectivity: 
=================================================================================
id: 4 epsilon: 5 mu: 3
vertices: 4 5 
subgraph: (cores: 4 5 non-cores: )
connectivity: (4,2):0.571429 
=================================================================================
id: 5 epsilon: 4 mu: 5
vertices: 4 
subgraph: (cores: 4 non-cores: 6 )
connectivity: (4,2):0.571429 
=================================================================================
id: 6 epsilon: 4 mu: 4
vertices: 5 8 
subgraph: (cores: 8 non-cores: 7 )(cores: 5 non-cores: )
connectivity: (5,2):0.571429 (8,10):0.6 
=================================================================================
id: 7 epsilon: 3 mu: 5
vertices: 5 
subgraph: (cores: 5 non-cores: 7 )
connectivity: (5,2):0.571429 
=================================================================================
id: 8 epsilon: 4 mu: 2
vertices: 6 7 
subgraph: (cores: 7 non-cores: )(cores: 6 non-cores: )
connectivity: (6,4):0.4 (7,8):0.4 
=================================================================================
id: 9 epsilon: 3 mu: 3
vertices: 7 
subgraph: (cores: 7 non-cores: )
connectivity: (7,8):0.4 (7,5):0.333333 
=================================================================================
id: 10 epsilon: 9 mu: 2
vertices: 9 10 
subgraph: (cores: 9 10 non-cores: )
connectivity: 
=================================================================================
id: 11 epsilon: 7 mu: 3
vertices: 9 10 11 
subgraph: (cores: 9 10 11 non-cores: )
connectivity: 
