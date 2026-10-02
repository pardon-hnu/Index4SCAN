# Two space-efficient indexes for high-speed structural graph clustering

> this repository is used for the paper "Advancing Structural Graph Clustering: Space-Efficient Indexes for High-Speed Queries"

## file structure

```
---datasets: consits of two toy graph, hiv and example. Example is the example graph in paper.
---src: source code
------example: the result of example, used in paper
------HeadFile: all head files
---------global: head files of global definitions and SCAN 
---------index: four indexes 
---------util: utils
------Index: storage of all indexes
------res: experiment results, used in paper
------Result: clustering results
---------lables: written by labels (binary file)
---------result: written by cluster sets
------shell_files: shell files for experiment in paper
------main.cpp
------makefile
------other auxiliary files
---VDStar: the source code provided by author of VD-STAR (KDD)
---VDStarNoT: the source code provided by author of VD-STAR (KDD)
```

## compile

```shell
cd src
make
```

## run

```shell
./main [Dataset] [Method] [Pattern] [Parameter 1]...[Parameter n]

Dataset:

Method:
    GS-Index: the sota index in query efficiency
    BOBTIN: the sota index in maintenance and construction
            rho: rho-approximate SCAN
            failure_pb: failure probability of bottom-k-sketch
            delta: the discretized number of similarity interval
    FOREST: forest index 
            delta: the discretized number of similarity interval
    PPT: ppt index 
            delta: the discretized number of similarity interval
            format_id:  the data structure of organizing PPT
                -1: SortSet
                -2: Table
                -3: RStar
Pattern:
    construct: index construction and storage
    query: a single query, following parameters: [epsilon] [mu]
    exp-query: query under varing parameters
    test: the pattern for testing 
    update: update, following parameters: ["insert" or "remove"] [u] [v] [0: rebuild 1: local reconstruction]
    exp-update: update under 1024 times ["insert" or "remove"] [0: rebuild 1: local reconstruction]
```

## Example:

```shell
    cd src
    make
    ./main example GS-Index query 0.6 5
    ./main example BOTBIN query 0.6 5
    ./main example FOREST query 0.6 5
    ./main example PPT query 0.6 5
    ./main hiv GS-Index query 0.6 5
    ./main hiv BOTBIN query 0.6 5
    ./main hiv FOREST query 0.6 5
    ./main hiv PPT query 0.6 5
```