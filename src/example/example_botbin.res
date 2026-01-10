[LOAD GRAPH] start!!! graph: example filepath: ../datasets/example.txt
graph size: 12
edge size: 18
[LOAD GRAPH] success!!!
[CONSTRUCT INDEX] start building BOTTOM-K-SKETCH
[CONSTRUCT INDEX] finish building BOTTOM-K-SKETCH
[TIME COST] 4.0379e-05 s.
[CONSTRUCT INDEX] start building Neighbor Order
[CONSTRUCT INDEX] finish building Neighbor Order
[TIME COST] 1.4033e-05 s.
[CONSTRUCT INDEX] start building bucket index
[CONSTRUCT INDEX] finish building bucket index
[TIME COST] 1.7936e-05 s.
[PRINT TIME] total building time: 7.2348e-05
[PRINT SPAEC]: total space cost: 0.00115204 MB
[PRINT HASH FOR EACH VERTEX] start !!!
0: 4
1: 2
2: 10
3: 9
4: 3
5: 7
6: 5
7: 11
8: 8
9: 6
10: 0
11: 1
[PRINT HASH FOR EACH VERTEX] finish !!!
[PRINT SKETCH] start !!!
0: 2 4 9 10 
1: 2 3 4 10 
2: 2 3 4 7 9 10 
3: 4 7 9 10 
4: 2 3 5 7 10 
5: 3 7 9 10 11 
6: 3 5 
7: 7 8 11 
8: 0 6 8 11 
9: 0 1 6 8 
10: 0 1 6 8 
11: 0 1 6 
[PRINT SKETCH] finish !!!
[PRINT NEIGHBOR ORDER] start !!!
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
[PRINT NEIGHBOR ORDER] finish !!!
[PRINT CLUSTER INDEX] start !!!
simialrity bucket [0,0.1): (5,2)(4,4)(4,5)(3,0)(3,1)(3,3)(3,8)(3,9)(3,10)(2,7)(2,11)(1,6)
simialrity bucket [0.1,0.2): (5,2)(4,4)(4,5)(3,0)(3,1)(3,3)(3,8)(3,9)(3,10)(2,7)(2,11)(1,6)
simialrity bucket [0.2,0.3): (5,2)(4,4)(4,5)(3,0)(3,1)(3,3)(3,8)(3,9)(3,10)(2,7)(2,11)(1,6)
simialrity bucket [0.3,0.4): (5,2)(4,4)(4,5)(3,0)(3,1)(3,3)(3,8)(3,9)(3,10)(2,7)(2,11)(1,6)
simialrity bucket [0.4,0.5): (5,2)(4,4)(3,0)(3,1)(3,3)(3,5)(3,8)(3,9)(3,10)(2,11)(1,6)(1,7)
simialrity bucket [0.5,0.6): (5,2)(3,0)(3,1)(3,3)(3,9)(3,10)(2,4)(2,5)(2,8)(2,11)
simialrity bucket [0.6,0.7): (3,0)(3,2)(3,9)(3,10)(2,1)(2,3)(2,8)(2,11)
simialrity bucket [0.7,0.8): (2,9)(2,10)(2,11)
simialrity bucket [0.8,0.9): (1,9)(1,10)
simialrity bucket [0.9,1): (1,9)(1,10)
[PRINT CLUSTER INDEX] finish !!!
