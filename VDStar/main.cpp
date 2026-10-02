#include "MyLib/MyTimer.h"
#include "MyLib/ParaReader.h"
#include "graph/Graph.h"
#include <iostream>
#include <set>
#include <stdio.h>
#include <random>
#include <fstream>

int n, m;
int **adj;
int *degree;

struct Param {
    char graph_file[200];
    double rho;
};

struct Param parseArgs(int nargs, char **args) {
    Param rtn;
    int cnt = 1;
    bool failed = false;
    char *arg;
    int i;
    char para[10];
    char graph_file[200];
    double rho = 0.01;

    printf("The input parameters are:\n\n");
    while (cnt < nargs && !failed) {
        arg = args[cnt++];
        if (cnt == nargs) {
            failed = true;
            break;
        }
        i = getNextChar(arg);
        if (arg[i] != '-') {
            failed = true;
            break;
        }
        getNextWord(arg + i + 1, para);
        printf("%s\t", para);
        arg = args[cnt++];
        if (strcmp(para, "graph") == 0) {
            getNextWord(arg, graph_file);
            printf("%s\n", graph_file);
        } else if (strcmp(para, "rho") == 0) {
            rho = atof(arg);
            if (rho < 0 || rho > 1) {
                failed = true;
                break;
            }
            printf("rho : %lf\n", rho);
        } 
        else {
            failed = true;
            printf("Unknown option -%s!\n\n", para);
        }
    }

    /*****************************************************************************/
    strcpy(rtn.graph_file, graph_file);
    rtn.rho = rho;
    return rtn;
}

void usage() {
    printf("Usage:\n");
    printf("dynstrclu -graph [graph file] -update [update file] "
           "-rho [\\rho]\n");
}

int generateRandomInt(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);

    return dis(gen);
}

void read_graph(string filepath)
{
    cout << "[LOAD GRAPH] start!!! filepath: " << filepath << endl;

    ifstream infile(filepath);
    string line;

    getline(infile, line);

    char *temp;
    temp = strtok(const_cast<char *>(line.c_str()), " ");
    n = atoi(temp);

    temp = strtok(NULL, " ");
    m = atoi(temp);

    cout << "graph size: " << n << endl;
    cout << "edge size: " << m << endl;

    adj = new int *[n];
    degree = new int[n];

    int count = 0;

    while (getline(infile, line))
    {
        if (line.size() == 0)
            break;

        temp = strtok(const_cast<char *>(line.c_str()), " ");
        int d = atoi(temp);

        degree[count] = d;
        adj[count] = new int[d];

        for (int i = 0; i < d; i++)
        {
            temp = strtok(NULL, " ");
            int v = atoi(temp);
            adj[count][i] = v;
        }

        count++;
    }

    infile.close();

    cout << "[LOAD GRAPH] success!!!" << endl;
}

int main(int argc, char **argv) {
    printf("Start to parse the arguments\n");
    Param para = parseArgs(argc, argv);
    printf("Arguments parsed.\n");

    read_graph(para.graph_file);

    MyVector<dynscan::Vertex *> _vList;
    _vList.reserve(n);

    for (int i = 0; i < n; i++) {
        dynscan::Vertex *newVertex = new dynscan::Vertex(i + 1);
        _vList.push_back(newVertex);
    }
    Graph graph(_vList, para.rho);
    printf("Graph generation finished with %d vertices and %d edges.\n", n,
           m );

    _vList.release_space();

    cout<<"Start to insert edges and construct index!!!"<<endl;
    double start = getCurrentTime();

    for (int u = 0; u < n; u++)
    {
        for (int j = 0; j < degree[u]; j++)
        {
            int v = adj[u][j];
            // 防止重复边（无向图）
            if (u < v)
            {
                if(u>n||v>n) {cout<<"Error: vertex id out of range! u: "<<u<<" v: "<<v<<endl; exit(1);}
                graph.insertEdge(u+1,v+1);
            }
        }
    }
    graph.create_table();
    printf("table created\n");
    double end = getCurrentTime();
    printf("Total time used after inserting construct index: %.9lf \n", end - start);
    
    for(int i = 0; i < n; i++){
        delete[] adj[i];   
    }
    delete[] adj;          
    delete[] degree;       
    
    double q_time = 0;
    // vector<double> epsilons={0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
    vector<double> epsilons={0.8,0.7,0.6,0.5,0.4,0.3,0.2};
    // vector<double> epsilons={
    // 	// 0.90,
    // 	// 0.89,0.88,0.87,0.86,
    // 	// 0.85,
    // 	// 0.84,0.83,0.82,0.81,
    // 	0.80,
    // 	0.79,0.78,0.77,0.76,
    // 	0.75,
    // 	0.74,0.73,0.72,0.71,
    // 	0.70,
    // 	0.69,0.68,0.67,0.66,
    // 	0.65,
    // 	0.64,0.63,0.62,0.61,
    // 	0.60,
    // 	0.59,0.58,0.57,0.56,
    // 	0.55,
    // 	0.54,0.53,0.52,0.51,
    // 	0.50,
    // 	0.49,0.48,0.47,0.46,
    // 	0.45,
    // 	0.44,0.43,0.42,0.41,
    // 	0.40,
    // 	0.39,0.38,0.37,0.36,
    // 	0.35,
    // 	0.34,0.33,0.32,0.31,
    // 	0.30,
    // 	0.29,0.28,0.27,0.26,
    // 	0.25,
    // 	0.24,0.23,0.22,0.21,
    // 	0.20,
    	// 0.19,0.18,0.17,0.16,
    	// 0.15,
    	// 0.14,0.13,0.12,0.11,
    	// 0.10
    // };
    // vector<double> epsilons = {
    //     0.8,
    //     0.7,
    //     0.6,
    //     0.5,
    //     0.4,
    //     0.3,
    //     0.2
    // };
    vector<int> mus={2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    for(auto eps:epsilons)
    {
        for(auto mu:mus)
        {
            cout<<"epsilon: "<<eps<<" mu: "<<mu<<endl;
            auto time_begin=getCurrentTime();
            q_time = graph.query(eps, mu);
            auto time_end=getCurrentTime();
            auto cluster_time=time_end-time_begin;
            cout<<"[TIME COST] "<<cluster_time<<" s."<<endl;
            cout<<"core time: "<<q_time<<" s."<<endl;
            cout<<"bfs time: "<<cluster_time-q_time<<" s."<<endl;
        }
    }
}
