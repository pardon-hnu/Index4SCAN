#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <fstream>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <queue>
#include <random>
#include <map>
#include <numeric>
#include <unordered_set>
#include <fstream>
using namespace std;
// Sample updates only when both endpoints have original degree <= this limit.
constexpr int max_update_degree = 5000;
map<string, string> dataset_2_SCAN_file = {
	//tiny graphs
	// {"example","../datasets/example.txt"},
	// {"hiv","../datasets/hiv.txt"},
	// {"dolphins","../datasets/dolphins.txt"},
	// {"contiguous","../datasets/contiguous.txt"},
	// {"club","../datasets/zachary_karate_club.txt"},
	// {"tribes","/home/hnu/Disk0/ParDon/graph/Tribes.txt"},
	//real graphs
	// {"amazon","/home/hnu/Disk0/ParDon/graph/amazon.txt"},
	// {"dblp","/home/hnu/Disk0/ParDon/graph/dblp.txt"},
	// {"youtube","/home/hnu/Disk0/ParDon/graph/youtube.txt"},
	// {"skitter","/home/hnu/Disk0/ParDon/graph/skitter.txt"},
	// {"pokec","/home/hnu/Disk0/ParDon/graph/pokec.txt"},
	// {"topcats","/home/hnu/Disk0/ParDon/graph/topcats.txt"},
	// {"livejournal1","/home/hnu/Disk0/ParDon/graph/livejournal1.txt"},
	// {"livejournal2","/home/hnu/Disk0/ParDon/graph/livejournal2.txt"},
	// {"orkut1","/home/hnu/Disk0/ParDon/graph/orkut1.txt"},
	// {"orkut2","/home/hnu/Disk0/ParDon/graph/orkut2.txt"},
	// {"indochina","/home/hnu/Disk0/ParDon/graph/indochina.txt"},
	// {"uk2002","/home/hnu/Disk0/ParDon/graph/uk2002.txt"},
	//real graphs with billion-scale directed edges
	// {"webbase","/home/hnu/Disk0/ParDon/graph/webbase.txt"},
	{"uk2005","/home/hnu/Disk0/ParDon/graph/uk2005.txt"},
	// real graphs with billon-scale undirected edges
	{"it2004","/home/hnu/Disk0/ParDon/graph/it2004.txt"},
	{"friendster","/home/hnu/Disk0/ParDon/graph/friendster.txt"},
	// {"sk2005","/home/hnu/Disk0/ParDon/graph/sk2005.txt"}
};
bool cmp_adj(int x, int y)
{
	return x>y;
}
class update_edges_generation
{
    private:
        int** graph;
        int* degree;
        int n = 0;
        int m = 0;
        int dmax = 0;
        string dataset;
        
        inline bool has_edge(int u, int v)
        {
            return binary_search(graph[u],graph[u] + degree[u],v,cmp_adj);
        }

        inline uint64_t edge_key(int u, int v)
        {
            if (u > v) swap(u, v);

            return (static_cast<uint64_t>(static_cast<uint32_t>(u)) << 32)
                | static_cast<uint32_t>(v);
        }

    public:
        update_edges_generation(string dataset)
        {   
            this->dataset = dataset;
            string filepath = dataset_2_SCAN_file[dataset];
            read_graph(filepath);
            sort_adj();
            for(int i=0;i<n;i++)
            {
                if(degree[i]>dmax) dmax=degree[i];
            }
        }
        ~update_edges_generation()
        {
            if (graph != nullptr)
            {
                for (int i = 0; i < n; ++i)
                {
                    delete[] graph[i];
                }

                delete[] graph;
                graph = nullptr;
            }

            delete[] degree;
            degree = nullptr;
        }
        void read_graph(string filepath)
        {
            cout << "[LOAD GRAPH] start!!! graph: "<<dataset<<" filepath: "<<filepath << endl;
            ifstream infile;
            infile.open(filepath);
            string line;
            getline(infile, line);
            char* temp;
            temp = strtok(const_cast<char*>(line.c_str()), " ");
            this->n = atoi(temp);
            cout<<"graph size: "<<n<<endl;
            temp = strtok(NULL, " ");
            this->m = atoi(temp);
            cout<<"edge size: "<<m<<endl;
            graph=new int*[n];
            degree=new int[n];
            int count=0;
            while (getline(infile, line))
            {
                if (line.size() == 0)
                {
                    break;
                }
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int d = atoi(temp);
                degree[count]=d;
                graph[count]=new int[d];
                for(int i=0;i<d;i++)
                {
                    temp = strtok(NULL, " ");
                    int v= atoi(temp);
                    graph[count][i]=v;
                }
                count++;
            }
            infile.close();
            cout << "[LOAD GRAPH] success!!!" << endl;
        }
        void sort_adj()
        {
            for(int i=0;i<n;i++)
            {
                sort(graph[i],graph[i]+degree[i],cmp_adj);
            }
        }
        vector<pair<int,int>> generate_insert_edges()
        {
            const int UPDATE_NUM = 1024;

            vector<pair<int,int>> update_edges;
            update_edges.reserve(UPDATE_NUM);

            unordered_set<uint64_t> selected;
            selected.reserve(UPDATE_NUM * 2);

            mt19937 rng(2026);
            uniform_int_distribution<int> vertex_dist(0, n - 1);

            while ((int)update_edges.size() < UPDATE_NUM)
            {
                int u = vertex_dist(rng);
                int v = vertex_dist(rng);

                if (u == v)
                    continue;

                if (degree[u] > max_update_degree || degree[v] > max_update_degree)
                    continue;

                if (has_edge(u, v))
                    continue;

                int a = min(u, v);
                int b = max(u, v);

                uint64_t key = edge_key(a, b);

                if (!selected.insert(key).second)
                    continue;

                update_edges.emplace_back(a, b);
            }

            return update_edges;
        }
        vector<pair<int,int>> generate_remove_edges()
        {
            const int UPDATE_NUM = 1024;

            vector<pair<int,int>> update_edges;
            update_edges.reserve(UPDATE_NUM);

            /*
            * prefix[u] 表示：
            * graph[0] ~ graph[u-1]
            * 一共有多少个邻接点。
            * 对无向图来说 prefix[n] = 2m
            */
            vector<long long> prefix(n + 1, 0);

            for (int u = 0; u < n; ++u)
            {
                prefix[u + 1] = prefix[u] + degree[u];
            }

            long long total_adj = prefix[n];

            cout << "total adjacency entries: "<< total_adj << endl;

            mt19937_64 rng(2026);

            uniform_int_distribution<long long> edge_dist( 0, total_adj - 1);

            unordered_set<uint64_t> selected;
            selected.reserve(UPDATE_NUM * 2);

            while ((int)update_edges.size() < UPDATE_NUM)
            {
                // 从所有邻接表位置中随机选一个
                long long pos = edge_dist(rng);
                /*
                * 找到：prefix[u] <= pos < prefix[u+1]
                */
                int u = upper_bound( prefix.begin(), prefix.end(), pos) - prefix.begin() - 1;
                int index =static_cast<int>(pos - prefix[u]);
                int v = graph[u][index];

                if (degree[u] > max_update_degree || degree[v] > max_update_degree)
                    continue;

                // 无向边统一成小编号在前
                int a = min(u, v);
                int b = max(u, v);
                uint64_t key = edge_key(a, b);

                // 避免同一条边被删除两次
                if (!selected.insert(key).second)
                    continue;

                update_edges.emplace_back(a, b);
            }

            return update_edges;
        }

        // ============================================================
        // 1024
        // u v
        // ...
        // ============================================================
        void write_update_edges(string filepath, const vector<pair<int,int>>& update_edges)
        {
            ofstream fout(filepath);
            if (!fout)
            {
                cerr << "cannot open output file: " << filepath << endl;
                exit(1);
            }
            fout << update_edges.size() << '\n';
            for (auto &e : update_edges)
            {
                fout << e.first << " " << e.second << '\n';
            }
            fout.close();
            cout << "[WRITE UPDATE EDGES] " << filepath << " success, num = " << update_edges.size()<< endl;
        }
};
int main()
{
    for(auto item: dataset_2_SCAN_file)
    {
        string dataset_name=item.first;
        cout<<"================================================================"<<endl;
        cout<<"Processing dataset: "<<dataset_name<<endl;
        update_edges_generation worker=update_edges_generation(dataset_name);
        vector<pair<int,int>> insert_edges=worker.generate_insert_edges();
        vector<pair<int,int>> remove_edges=worker.generate_remove_edges();
        worker.write_update_edges("./updates/"+dataset_name+"/insert_1024.txt", insert_edges);
        worker.write_update_edges("./updates/"+dataset_name+"/remove_1024.txt", remove_edges);
    }
    return 0;
}
