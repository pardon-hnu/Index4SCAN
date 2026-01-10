[LOAD GRAPH] start!!! graph: example filepath: ../datasets/example.txt
graph size: 12
edge size: 18
[LOAD GRAPH] success!!!
[CONSTRUCT INDEX] start building BOTTOM-K-SKETCH
[CONSTRUCT INDEX] finish building BOTTOM-K-SKETCH
[TIME COST] 3.2829e-05 s.
[CONSTRUCT INDEX] start building Neighbor Order
[CONSTRUCT INDEX] finish building Neighbor Order
[TIME COST] 1.1401e-05 s.
[CONSTRUCT FOREST] start
[CONSTRUCT FOREST] now construct 0 th tree
[CONSTRUCT FOREST] now construct 1 th tree
[CONSTRUCT FOREST] now construct 2 th tree
[CONSTRUCT FOREST] now construct 3 th tree
[CONSTRUCT FOREST] now construct 4 th tree
[CONSTRUCT FOREST] now construct 5 th tree
[CONSTRUCT FOREST] now construct 6 th tree
[CONSTRUCT FOREST] now construct 7 th tree
[CONSTRUCT FOREST] now construct 8 th tree
[CONSTRUCT FOREST] finishi! time cost: 0.000354214
[PRINT NODE NUM] total node num of forest is 46
[PRINT TIME] total building time: 0.000398444
[PRINT SPACE INFORMATION] core vertices space cost: 0.000278473.MB
[PRINT SPACE INFORMATION] non core vertices space cost: 0.000221252.MB
[PRINT SPACE]: total space cost: 0.0011673 MB
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
[PRINT FOREST] start!!!
print the 0 th tree corresponding range (1,0.9]
node ID: 0 mu: 1 father: 0 childs: 1  cores: 
non cores: 
node ID: 1 mu: 2 father: 0 childs:  cores: 10 9 
non cores: 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 1, mu: 2  cores: 10 9   non cores: 
print the 1 th tree corresponding range (0.9,0.8]
node ID: 0 mu: 1 father: 0 childs: 1  cores: 
non cores: 
node ID: 1 mu: 2 father: 0 childs:  cores: 10 9 
non cores: 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 1, mu: 2  cores: 10 9   non cores: 
print the 2 th tree corresponding range (0.8,0.7]
node ID: 0 mu: 1 father: 0 childs: 1  cores: 
non cores: 
node ID: 1 mu: 3 father: 0 childs:  cores: 11 10 9 
non cores: 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 1, mu: 3  cores: 11 10 9   non cores: 
print the 3 th tree corresponding range (0.7,0.6]
node ID: 0 mu: 1 father: 0 childs: 3 4  cores: 
non cores: 
node ID: 4 mu: 3 father: 0 childs: 2  cores: 11 8 
non cores: 
node ID: 3 mu: 3 father: 0 childs: 1  cores: 3 1 
non cores: 
node ID: 2 mu: 4 father: 4 childs:  cores: 10 9 
non cores: 11 8 
node ID: 1 mu: 4 father: 3 childs:  cores: 2 0 
non cores: 3 1 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 3, mu: 3  cores: 3 1   non cores: 
    Node ID: 1, mu: 4  cores: 2 0   non cores: 3 1 
  Node ID: 4, mu: 3  cores: 11 8   non cores: 
    Node ID: 2, mu: 4  cores: 10 9   non cores: 11 8 
print the 4 th tree corresponding range (0.6,0.5]
node ID: 0 mu: 1 father: 0 childs: 4 5  cores: 
non cores: 
node ID: 5 mu: 3 father: 0 childs: 3  cores: 11 8 
non cores: 
node ID: 4 mu: 3 father: 0 childs: 2  cores: 5 4 
non cores: 
node ID: 3 mu: 4 father: 5 childs:  cores: 10 9 
non cores: 11 8 
node ID: 2 mu: 4 father: 4 childs: 1  cores: 3 1 0 
non cores: 5 4 
node ID: 1 mu: 6 father: 2 childs:  cores: 2 
non cores: 5 4 3 1 0 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 4, mu: 3  cores: 5 4   non cores: 
    Node ID: 2, mu: 4  cores: 3 1 0   non cores: 5 4 
      Node ID: 1, mu: 6  cores: 2   non cores: 5 4 3 1 0 
  Node ID: 5, mu: 3  cores: 11 8   non cores: 
    Node ID: 3, mu: 4  cores: 10 9   non cores: 11 8 
print the 5 th tree corresponding range (0.5,0.4]
node ID: 0 mu: 1 father: 0 childs: 6 7  cores: 
non cores: 
node ID: 7 mu: 2 father: 0 childs: 5  cores: 7 
non cores: 
node ID: 6 mu: 2 father: 0 childs: 3  cores: 6 
non cores: 
node ID: 5 mu: 3 father: 7 childs: 4  cores: 11 
non cores: 
node ID: 4 mu: 4 father: 5 childs:  cores: 10 9 8 
non cores: 11 7 
node ID: 3 mu: 4 father: 6 childs: 2  cores: 5 3 1 0 
non cores: 
node ID: 2 mu: 5 father: 3 childs: 1  cores: 4 
non cores: 6 5 1 
node ID: 1 mu: 6 father: 2 childs:  cores: 2 
non cores: 5 4 3 1 0 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 6, mu: 2  cores: 6   non cores: 
    Node ID: 3, mu: 4  cores: 5 3 1 0   non cores: 
      Node ID: 2, mu: 5  cores: 4   non cores: 6 5 1 
        Node ID: 1, mu: 6  cores: 2   non cores: 5 4 3 1 0 
  Node ID: 7, mu: 2  cores: 7   non cores: 
    Node ID: 5, mu: 3  cores: 11   non cores: 
      Node ID: 4, mu: 4  cores: 10 9 8   non cores: 11 7 
print the 6 th tree corresponding range (0.4,0.3]
node ID: 0 mu: 1 father: 0 childs: 6  cores: 
non cores: 
node ID: 6 mu: 2 father: 0 childs: 5  cores: 6 
non cores: 
node ID: 5 mu: 3 father: 6 childs: 3 4  cores: 11 7 
non cores: 
node ID: 4 mu: 4 father: 5 childs:  cores: 10 9 8 
non cores: 11 7 
node ID: 3 mu: 4 father: 5 childs: 2  cores: 3 1 0 
non cores: 
node ID: 2 mu: 5 father: 3 childs: 1  cores: 5 4 
non cores: 7 6 3 1 
node ID: 1 mu: 6 father: 2 childs:  cores: 2 
non cores: 5 4 3 1 0 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 6, mu: 2  cores: 6   non cores: 
    Node ID: 5, mu: 3  cores: 11 7   non cores: 
      Node ID: 3, mu: 4  cores: 3 1 0   non cores: 
        Node ID: 2, mu: 5  cores: 5 4   non cores: 7 6 3 1 
          Node ID: 1, mu: 6  cores: 2   non cores: 5 4 3 1 0 
      Node ID: 4, mu: 4  cores: 10 9 8   non cores: 11 7 
print the 7 th tree corresponding range (0.3,0.2]
node ID: 0 mu: 1 father: 0 childs: 6  cores: 
non cores: 
node ID: 6 mu: 2 father: 0 childs: 5  cores: 6 
non cores: 
node ID: 5 mu: 3 father: 6 childs: 3 4  cores: 11 7 
non cores: 
node ID: 4 mu: 4 father: 5 childs:  cores: 10 8 
non cores: 11 7 
node ID: 3 mu: 4 father: 5 childs: 2  cores: 9 3 1 0 
non cores: 11 
node ID: 2 mu: 5 father: 3 childs: 1  cores: 5 4 
non cores: 7 6 3 1 
node ID: 1 mu: 6 father: 2 childs:  cores: 2 
non cores: 5 4 3 1 0 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 6, mu: 2  cores: 6   non cores: 
    Node ID: 5, mu: 3  cores: 11 7   non cores: 
      Node ID: 3, mu: 4  cores: 9 3 1 0   non cores: 11 
        Node ID: 2, mu: 5  cores: 5 4   non cores: 7 6 3 1 
          Node ID: 1, mu: 6  cores: 2   non cores: 5 4 3 1 0 
      Node ID: 4, mu: 4  cores: 10 8   non cores: 11 7 
print the 8 th tree corresponding range (0.2,0.1]
node ID: 0 mu: 1 father: 0 childs: 6  cores: 
non cores: 
node ID: 6 mu: 2 father: 0 childs: 5  cores: 6 
non cores: 
node ID: 5 mu: 3 father: 6 childs: 3 4  cores: 11 7 
non cores: 
node ID: 4 mu: 4 father: 5 childs:  cores: 10 8 
non cores: 11 7 
node ID: 3 mu: 4 father: 5 childs: 2  cores: 9 3 1 0 
non cores: 11 
node ID: 2 mu: 5 father: 3 childs: 1  cores: 5 4 
non cores: 7 6 3 1 
node ID: 1 mu: 6 father: 2 childs:  cores: 2 
non cores: 5 4 3 1 0 
Node ID: 0, mu: 1  cores:   non cores: 
  Node ID: 6, mu: 2  cores: 6   non cores: 
    Node ID: 5, mu: 3  cores: 11 7   non cores: 
      Node ID: 3, mu: 4  cores: 9 3 1 0   non cores: 11 
        Node ID: 2, mu: 5  cores: 5 4   non cores: 7 6 3 1 
          Node ID: 1, mu: 6  cores: 2   non cores: 5 4 3 1 0 
      Node ID: 4, mu: 4  cores: 10 8   non cores: 11 7 
print the 9 th tree corresponding range (0.1,0]
[PRINT FOREST] finish!!!
