# Two space-efficient indexes for high-speed structural graph clustering
> this respository is used for the paper "Advancing Structural Graph Clustering: Space-Efficient Indexes for High-Speed Queries"

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
                -3: PST
                -4: KD
                -5: RStar
Pattern:
    construct: index construction and storage
    query: a single query, following parameters: [epsilon] [mu]
    exp-query: query under varing parameters
    test: the pattern for testing 
```
## Example:
```shell
    cd src
    make
    ./main hiv GS-Index construct
    ./main hiv GS-Index query 0.6 5
    ./main hiv BOTBIN construct
    ./main hiv BOTBIN query 0.6 5
    ./main hiv FOREST construct
    ./main hiv FOREST query 0.6 5
    ./main hiv PPT construct
    ./main hiv PPT query 0.6 5
```