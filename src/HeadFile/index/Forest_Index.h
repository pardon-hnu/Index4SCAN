#pragma once

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
#include <numeric>
#include <unordered_set>
#include <climits>
#include <global/global.h>
#include <index/Index.h>
#include <util/RandList.h>
#include <util/SDSU.h>
using namespace std;
// #define DEBUG 
// #define PRINT 
// #define DEBUG_CLUSTERING
bool cmp_order_forest(pair<float,int> x,pair<float,int> y)
{
	return x.first>y.first;
}
bool cmp_order_priority(pair<int,int> x,pair<int,int> y)
{
	return x.second>y.second;
}
class FNode    
{
    public:
        int mu;
        int ID;
        int father;
        vector<int> cores;
        vector<int> non_cores;
        vector<int> childs;
        FNode()
        {

        }
        FNode(int mu,int ID)
        {
            this->mu=mu;
            this->ID=ID;
            this->father=ID;
        }
};

typedef vector<FNode> tree_t;

typedef struct
{
    int u;
    int v;
    int old_similarity;
    int new_similarity;
    int update_type;
}update_edge_in_graph_t;
typedef struct
{
    int u;
    int v;
    int old_level;
    int new_level;
}update_edge_in_tree_t;

typedef struct
{
    vector<pair<int,int>> new_core_v;
    vector<int> new_connectiviy_v;
    vector<update_edge_in_graph_t> update_edges;
}tree_update_info_t;

class SCAN_Forest_Index : public SCAN_Index
{
    private:
        int k=0;
        float RHO=0.1;
        float FAILURE_PB=0.001;
        int DELTA=10;
        static vector<int> rnk;
        vector<vector<int> > sketch;

        SCAN_DSU dsu;
        UnionFind update_uf;

        vector<vector<pair<int, int>>> neighbor_order;
        vector<tree_t> forest;

        vector<vector<int>> vertex2node;

        int forest_update_token=0;
        vector<int> affected_vertex_mark;
        vector<int> active_vertex_mark;
        vector<int> ct_delta_mark;
        vector<int> update_ct_delta;
        vector<int> affected_node_mark;
        vector<int> affected_root_mark;
        vector<int> old_root_cache_mark;
        vector<int> old_root_cache;

        int total_node_num=0;
        int traverse_node_num=0;
        int traverse_core_num=0;
        int traverse_non_core_num=0;

        double core_vertices_space=0.0;
        double non_core_vertices_space=0.0;

        double total_ari=0.0;
        int query_count=0;
        int update_count=0;
    public:
        SCAN_Forest_Index(string dataset,int DELTA=10)
        {
            this->DELTA=DELTA;
            forest.resize(DELTA);
            init_graph(dataset);
            k=(int)(ceil((1/(2*pow(RHO,2)))*(log(2/FAILURE_PB))));
            dsu=SCAN_DSU(n);
            update_uf=UnionFind(n);
            affected_vertex_mark.resize(n,0);
            active_vertex_mark.resize(n,0);
            ct_delta_mark.resize(n,0);
            update_ct_delta.resize(n,0);
        }
        ~SCAN_Forest_Index()
        {

        }
        static bool cmp_rnk(const int& u, const int& v) {
            return rnk[u] < rnk[v];
        }
        int find_kthelement_ab(const vector<int>& a, const vector<int>& b, const int& k) {
            if (k > a.size() + b.size()) 
                return INF;
            int cnt = 0, res = 0;
            auto cur1 = a.begin();
            auto cur2 = b.begin();
            while (cnt != k) {
                if (*cur1 < *cur2) {
                    cnt++;
                    res = *cur1;
                    cur1++;
        
                }
                else if (*cur1 == *cur2) {
                    cnt++;
                    res = *cur1;
                    cur1++;
                    cur2++;
                }
                else {
                    cnt++;
                    res = *cur2;
                    cur2++;
                }
                if (cur1 == a.end() || cur2 == b.end())
                    break;
            }
            //cout << cnt << "\n";
        
            if (cnt < k) {
                if (cur1 == a.end()) {
                    for (; cnt != k && cur2 != b.end(); cnt++) {
                        res = *cur2;
                        cur2++;
                    }
                }
                else {
                    for (; cnt != k && cur1 != a.end(); cnt++) {
                        res = *cur1;
                        cur1++;
                    }
                }
            }
            if (cnt < k)
                return INF;
            return res;
        }
        int find_common_ab(vector<int> & a, vector<int> & b)
        {
            auto cur_a=a.begin();
            auto cur_b=b.begin();
            int count=0;
            while(true)
            {
                if (cur_a == a.end() || cur_b == b.end())
                    break;
                if(*cur_a==*cur_b)
                {
                    count++;
                    cur_a++;
                }
                else if(*cur_a>*cur_b)
                {
                    *cur_b++;
                }
                else
                {
                    *cur_a++;
                }
            }
            return count;
        }
        int find_common_ab_k(vector<int> & a, vector<int> & b,int cur_k)
        {
            auto cur_a=a.begin();
            auto cur_b=b.begin();
            int count=0;
            while(true)
            {
                if (cur_a == a.end() || cur_b == b.end())
                    break;
                if(*cur_a==*cur_b)
                {
                    
                    cur_a++;
                    if(*cur_a<=cur_k) count++;
                    else break;
                }
                else if(*cur_a>*cur_b)
                {
                    *cur_b++;
                }
                else
                {
                    *cur_a++;
                }
            }
            return count;
        }        
        int compute_approximate_jaccard_similarity(int u, int v)
        {
            float j_sim=0;
            int cn=0;
            int cur_k=find_kthelement_ab(sketch[u], sketch[v], k);
            
            if (cur_k == INF) {
				cn=find_common_ab(sketch[u],sketch[v]);
				j_sim = float(cn) / (sketch[u].size()+sketch[v].size()-cn);
			}
			else {
                cn=find_common_ab_k(sketch[u],sketch[v],cur_k);
				j_sim = float(cn) / k;
			}
            return j_sim*PRECISION;
        }
        void construct_bottom_k_sketch()
        {
            cout<<"[CONSTRUCT INDEX] start building BOTTOM-K-SKETCH"<<endl;
            const auto begin=std::chrono::steady_clock().now();

            mt19937 rng(random_device{}());
            rnk.resize(n);
            iota(rnk.begin(), rnk.end(), 0);
            shuffle(rnk.begin(), rnk.end(), rng);

            int* rnk_id = new int[n];
            for (int i = 0; i < n; i++) {
                rnk_id[i] = i;
            }
            sort(rnk_id, rnk_id + n, cmp_rnk);
            sketch.resize(n);

            for (int i = 0; i < n; i++) {

                int u = rnk_id[i];
                if (sketch[u].size() < k)
                    sketch[u].push_back(i);

                for (int j = 0; j < degree[u]; j++) {
                    int v = graph[u][j];
                    if (sketch[v].size() >= k) {
                        continue;
                    }
                    sketch[v].push_back(i);
                }
            }
            for(int i=0;i<n;i++)
            {
                sort(sketch[i].begin(),sketch[i].end());
            }

            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish building BOTTOM-K-SKETCH"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            total_build_time+=construct_time.count();
        }  
        float jaccard_similarity(int common, int lenth1,int lenth2)
        {
            float score= (float)(common)/(float)(lenth1+lenth2-common);
            return score;
        }
        int sort_and_merge_set_intersection_cardinality(int * set1,int lenth1, int * set2, int lenth2)
        {
            int i1=0;
            int i2=0;
            int count=0;
            while(true)
            {
                if(i1>=lenth1||i2>=lenth2)
                {
                    break;
                }
                if(set1[i1]==set2[i2])
                {
                    count++;
                    i1++;
                }
                else if(set1[i1]>set2[i2])
                {
                    i1++;
                }
                else
                {
                    i2++;
                }
            }
            return count;
        }
        int compute_similarity_exact(int v1,int v2)
        {
            int common=sort_and_merge_set_intersection_cardinality(graph[v1],degree[v1],graph[v2],degree[v2]);
            common+=2;
            float score=jaccard_similarity(common,degree[v1]+1,degree[v2]+1);
            return score*PRECISION;
        }
        void construct_neighbor_order()
        {
            cout<<"[CONSTRUCT INDEX] start building Neighbor Order"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            neighbor_order.resize(n);
            for(int u=0;u<n;u++)
            {
                for(int j=0;j<degree[u];j++)
                {
                    int v=graph[u][j];
                    if(v>u)
                    {
                        int similarity=compute_approximate_jaccard_similarity(u,v);
                        neighbor_order[u].push_back(make_pair(similarity,v));
                        neighbor_order[v].push_back(make_pair(similarity,u));
                    }
                }
                sort(neighbor_order[u].begin(),neighbor_order[u].end(),cmp_order_forest);
            }
            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish building Neighbor Order"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            total_build_time+=construct_time.count();
        }
        void fast_construct_neighbor_order()
        {
            string neighbor_order_file="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            ifstream read_index;
            read_index.open(neighbor_order_file);
            string line;
            while(getline(read_index, line))
            {
                char* temp;
                vector<pair<int,int>> item;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int len=atoi(temp);
                for(int i=0;i<len;i++)
                {
                    temp = strtok(NULL, " ");
                    int sim=atoi(temp);
                    temp = strtok(NULL, " ");
                    int v=atoi(temp);
                    item.push_back(make_pair(sim,v));
                }
                neighbor_order.push_back(item);
            }
            read_index.close();

            string construct_neighbor_order_time_file="./res/construct/"+dataset+"/botbin.txt";
            double construct_time = 0.0;
            bool read_time=true;

            read_index.open(construct_neighbor_order_time_file);
            while(getline(read_index, line))
            {
                if (line.find("start building BOTTOM-K-SKETCH") != string::npos || line.find("start building Neighbor Order") != string::npos) 
                {
                    read_time=true;
                }
                if (line.find("[TIME COST]") != string::npos && read_time) 
                {
                    // extract number
                    size_t pos1 = line.find("]");
                    size_t pos2 = line.find("s");
                    if (pos1 == string::npos || pos2 == string::npos) continue;
                    double t = atof(line.substr(pos1 + 1, pos2 - pos1 - 1).c_str());
                    // cout<<t<<endl;
                    construct_time += t;
                    read_time=false;
                }
            }
            read_index.close();
            cout<<"[CONSTRUCT INDEX] finish loading Neighbor Order"<<endl;
            cout<<"[TIME COST] "<<construct_time<<" s."<<endl;
            total_build_time+=construct_time;
        }
        void construct_forest()
        {
            cout<<"[CONSTRUCT FOREST] start"<<endl;
            auto construct_forest_begin=std::chrono::steady_clock().now();
            if(vertex2node.size()<DELTA) vertex2node.resize(DELTA);
            vector<int> subgraph;
            subgraph.reserve(n);
            vector<bool> vertex_bit_set;
            vertex_bit_set.resize(n);
            std::fill_n(vertex_bit_set.begin(), vertex_bit_set.size(), false);
            UnionFind UF(n);
            for(int tree_id=DELTA-1;tree_id>0;--tree_id)
            {
                cout<<"[CONSTRUCT FOREST] now construct "<<tree_id<<" th tree"<<endl;
                int range=tree_id*(PRECISION/DELTA);

                #ifdef DEBUG
                    cout<<"init subgraph, vertex_bit_set, tree, root_node, and verte2node map"<<endl;
                #endif
                //init 
                UF.Init(subgraph);
                for(auto v:subgraph)
                {vertex_bit_set[v]=false;}

                //init i-th tree and root node
                tree_t tree;
                FNode root_node=FNode(1,0);
                tree.push_back(root_node);
                int ID=1;

                //init a map: vertex id->tree node id
                vector<int> vertex_node_map(n,-1);

                #ifdef DEBUG
                    cout<<"core threshold decomposition"<<endl;
                #endif
                //get range-neighbors number for all vertices
                vector<vector<int>> v_set_of_each_core_threshold;
                v_set_of_each_core_threshold.resize(dmax+2);
                for(int i=0;i<n;i++)
                {
                    int core_threshold=1;
                    for(int j=0;j<neighbor_order[i].size();j++)
                    {
                        if(neighbor_order[i][j].first<range) break;
                        core_threshold++;
                    }
                    v_set_of_each_core_threshold[core_threshold].push_back(i);
                }

                //bottom up
                for(int temp_mu=dmax+1;temp_mu>1;--temp_mu)
                {
                    //ensure exists (epsilon,mu)-level
                    if(v_set_of_each_core_threshold[temp_mu].size()==0) continue;
                    #ifdef DEBUG
                    // cout<<"core threshold: "<<temp_mu<<endl;
                    #endif
                    /*
                        get the child tree node for each vertex v
                    */
                    unordered_map<int,vector<int>> vertex_child_node_map;
                    for(auto v:v_set_of_each_core_threshold[temp_mu])
                    {
                        vector<int> child_nodes;
                        for(auto &nv:neighbor_order[v])
                        {
                            //if nv is in the higher-level 
                            if(nv.first<range) break;
                            if(vertex_bit_set[nv.second])
                            {
                                int child_node_id=vertex_node_map[nv.second];
                                tree[child_node_id].non_cores.push_back(v);
                                int root_child_node_id=child_node_id;
                                while(tree[root_child_node_id].father!=0)
                                {
                                    root_child_node_id=tree[root_child_node_id].father;
                                }
                                child_nodes.push_back(root_child_node_id);
                            }
                        }
                        vertex_child_node_map[v]=child_nodes;
                    }

                    //add v[mu] into subgraph
                    for(auto v:v_set_of_each_core_threshold[temp_mu])
                    {
                        subgraph.push_back(v);
                        vertex_bit_set[v]=true;
                    }

                    //connect vertices
                    for(auto v:v_set_of_each_core_threshold[temp_mu])
                    {
                        for(auto &nv:neighbor_order[v])
                        {
                            if(nv.first<range) break;
                            if(vertex_bit_set[nv.second]&&UF.Find(v)!=UF.Find(nv.second))
                            {UF.Unite(v,nv.second);}
                        }
                    }

                    //divided v[mu] into different node
                    unordered_map<int,vector<int>> root2vertex;
                    for(auto v:v_set_of_each_core_threshold[temp_mu])
                    {
                        root2vertex[UF.Find(v)].push_back(v);
                    }
                    for(auto item:root2vertex)
                    {
                        FNode temp_node=FNode(temp_mu,ID);
                        temp_node.father=0;
                        for(auto v:item.second)
                        {
                            temp_node.cores.push_back(v);
                            vertex_node_map[v]=ID;
                            for(auto child:vertex_child_node_map[v])
                            {
                                temp_node.childs.push_back(child);
                                tree[child].father=ID;
                            }
                        }
                        tree.push_back(temp_node);
                        ID++;
                    }
                    
                }
                for(auto &item:tree)
                {
                    if(item.father==0&&item.ID!=0) tree[0].childs.push_back(item.ID);
                    sort(item.cores.begin(), item.cores.end());
                    item.cores.erase(unique(item.cores.begin(), item.cores.end()), item.cores.end());
                    sort(item.non_cores.begin(), item.non_cores.end());
                    item.non_cores.erase(unique(item.non_cores.begin(), item.non_cores.end()), item.non_cores.end());
                    sort(item.childs.begin(), item.childs.end());
                    item.childs.erase(unique(item.childs.begin(), item.childs.end()), item.childs.end());
                }
                forest[tree_id]=move(tree);
                vertex2node[tree_id]=move(vertex_node_map);
            }
            auto construct_forest_end=std::chrono::steady_clock().now();
            tms cosntruct_forest_time=construct_forest_end-construct_forest_begin;
            cout<<"[CONSTRUCT FOREST] finishi! time cost: "<<cosntruct_forest_time.count()<<endl;
            total_build_time+=cosntruct_forest_time.count();
            
            for(auto &tree:forest)
            {
                for(auto &item:tree)
                {
                    core_vertices_space+=item.cores.size()*4;
                    non_core_vertices_space+=item.non_cores.size()*4;
                    total_build_space+=(3+item.childs.size()+item.cores.size()+item.non_cores.size())*4;
                }
                total_node_num+=tree.size();
            }
        }  
        void construct() override
        {
            #ifdef FAST_CONSTRUCT
                fast_construct_neighbor_order();
            #else
                construct_bottom_k_sketch();
                construct_neighbor_order();
            #endif
            construct_forest();
        }
        int locate(int epsilon)
        {
            int bucket_range=PRECISION/DELTA;
            int bucket_id= epsilon/bucket_range;
            return bucket_id;
        }
        void traverse_tree_query(const tree_t&tree, int & mu, int nodeID, int & cluster_num)
        {
            if (nodeID >= tree.size()) return;
            // traverse_node_num++;
            const FNode & node=tree[nodeID];
            if(node.mu>=mu)
            {
                queue<int> q;
                q.push(node.ID);
                while(!q.empty())
                {
                    int temp_id=q.front();
                    q.pop();
                    const FNode & temp_node=tree[temp_id];
                    for(auto core:temp_node.cores)
                    {
                        labels[core]=cluster_num;
                        // traverse_core_num++;
                    }
                    for(auto non_core:temp_node.non_cores)
                    {
                        labels[non_core]=cluster_num;
                        // traverse_non_core_num++;
                    } 
                    for(auto child:temp_node.childs)
                    {q.push(child);}
                }
                cluster_num++;
            }
            else
            {
                for(auto child:node.childs)
                {
                    traverse_tree_query(tree,mu,child,cluster_num);
                }
            }
        }
        void forest_query(int epsilon, int mu)
        {
            clear_cluster_result();
            // traverse_node_num=0;
            // traverse_core_num=0;
            // traverse_non_core_num=0;
            cout<<"[QUERY] Start"<<endl;
            auto query_begin=std::chrono::steady_clock().now();
            int locate_tree=locate(epsilon);
            cout<<"search from "<<locate_tree<<" th tree"<<endl;
            int cluster_num=0;
            traverse_tree_query(forest[locate_tree],mu,0,cluster_num);
            auto query_end=std::chrono::steady_clock().now();
            tms query_time=query_end-query_begin;
            // cout<<"traverse node num: "<<traverse_node_num<<endl;
            // cout<<"traverse core num: "<<traverse_core_num<<endl;
            // cout<<"traverse non-core num: "<<traverse_non_core_num<<endl;
            cout<<"[QUERY] finish"<<endl;
            cout<<"[TIME COST] "<<query_time.count()<<" s."<<endl;
            cout<<"[QUERY] cluster number: "<<cluster_num<<endl;
            total_cluster_time+=query_time.count();
            query_count++;
        }
        void query(float epsilon, int mu) override
        {
            int eps=epsilon*PRECISION;
            #ifdef DEBUG_CLUSTERING
                check_query(eps,mu);
            #else
                forest_query(eps,mu);
            #endif
            #ifdef EVALUATION_CLUSTERING_QUALITY
                evaluate_clustering_quality(epsilon,mu);
            #endif
        }
        void evaluate_clustering_quality(float epsilon, int mu) 
        {
            // string groundtruth_file="Result/labels/"+dataset+"/GS-Index-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".txt";
            string groundtruth_file="Result/labels/"+dataset+"/GS-Index-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".bin";
            // vector<int> groundtruth;
            // groundtruth.reserve(n);
            // ifstream read_groundtruth;
            // read_groundtruth.open(groundtruth_file);
            // string line;
            // getline(read_groundtruth, line);
            // char* temp;
            // temp = strtok(const_cast<char*>(line.c_str()), " ");
            // while (temp != NULL) {
            //     int label = atoi(temp);
            //     groundtruth.push_back(label);
            //     temp = strtok(NULL, " ");
            // }
            // read_groundtruth.close();
            ifstream read_groundtruth(groundtruth_file, ios::binary);
            vector<int> groundtruth(n);
            read_groundtruth.read(reinterpret_cast<char*>(groundtruth.data()),sizeof(int)*n);
            read_groundtruth.close();
            double temp_ari=adjusted_rand_index_ignore_minus1(groundtruth,labels);
            cout<<"[EVALUATION] ari index score: "<<temp_ari<<endl;
            total_ari+=temp_ari;
        }
        void printTree(const tree_t& tree) 
        {
            for(auto &node:tree)
            {
                cout<<"node ID: "<<node.ID<<" ";
                cout<<"mu: "<<node.mu<<" ";
                cout<<"father: "<<node.father<<" ";
                cout<<"childs: ";
                for(auto &c:node.childs)
                {
                    cout<<c<<' ';
                }
                cout<<" ";
                cout<<"cores: ";
                for(auto &v:node.cores)
                {
                    cout<<v<<' ';
                }
                cout<<endl;
                cout<<"non cores: ";
                for(auto &nc:node.non_cores)
                {
                    cout<<nc<<' ';
                }
                cout<<endl;
            }
        }
        void traverse_print_tree(const tree_t& tree, int nodeID, int depth = 0)
        {
            if (nodeID >= tree.size()) return;

            const FNode& node = tree.at(nodeID);

            // 打印缩进
            string indent(depth * 2, ' ');
            cout << indent << "Node ID: " << node.ID << ", mu: " << node.mu;

            // 打印 cores
            cout<< "  cores: ";
            for (int v : node.cores) {
                cout << v << " ";
            }
            cout<<"  non cores: ";
            for(int v:node.non_cores)
            {
                cout<<v<<' ';
            }
            cout << endl;
            

            // 递归打印子节点
            if (node.childs.empty()) return;
            for (auto childID : node.childs) 
            {
                traverse_print_tree(tree, childID, depth + 1);
            }
        }
        void print_build_cost() override
        {
            cout<<"=================================================================="<<endl;
            cout<<"[PRINT NDDE NUM] each node num of forest is :"<<endl;
            for(int i=0;i<DELTA;i++)
            {
                cout<<i<<" th tree node num: "<<forest[i].size()<<endl;
            }
            cout<<"[PRINT NODE NUM] total node num of forest is "<<total_node_num<<endl;
            cout<<"[PRINT NDDE NUM] average node num of each tree is "<<double(total_node_num)/double(DELTA)<<endl;
            cout<<"[PRINT TIME] total building time: "<<total_build_time<<endl;
            cout<<"[PRINT SPACE INFORMATION] core vertices space cost: "<<(core_vertices_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] non core vertices space cost: "<<(non_core_vertices_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE]: total space cost: "<<(total_build_space/1024)/1024<<" MB"<<endl;
        }
        void print_cluster_time() override
        {
            cout<<"[PRINT TIME] total clustering time: "<<total_cluster_time<<endl;
            cout<<"[PRINT TIME] average clustering time: "<<total_cluster_time/double(query_count)<<endl;
            #ifdef EVALUATION_CLUSTERING_QUALITY
                cout<<"[PRINT ARI] average ARI score: "<<total_ari/double(query_count)<<endl;
            #endif
        }
        void print_index() override
        {
            cout<<"[PRINT NEIGHBOR ORDER] start !!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<i<<": ";
                for(auto &item:neighbor_order[i])
                {
                    cout<<"("<<item.first<<","<<item.second<<")";
                }
                cout<<endl;
            }
            cout<<"[PRINT NEIGHBOR ORDER] finish !!!"<<endl;
            cout<<"[PRINT FOREST] start!!!"<<endl;
            for(int i=0;i<DELTA;i++)
            {
                float num=DELTA;
                cout<<"print the "<<i<<" th tree corresponding range ["<<i/num<<","<<(i+1)/num<<") "<<endl;
                printTree(forest[i]);
                cout<<"print the struct of tree"<<endl;
                traverse_print_tree(forest[i],0);
            }
            // int i=7;
            // cout<<"print the "<<i<<" th tree corresponding range ("<<1.00-float(i)/float(DELTA)<<","<<1.00-float(i+1)/float(DELTA)<<"]"<<endl;
            // printTree(forest[i]);
            // traverse_print_tree(forest[i],0);
            cout<<"[PRINT FOREST] finish!!!"<<endl;
        }
        void store_index() override
        {
            cout<<"[STORE INDEX] start !!!!"<<endl;
            string store_forest="./Index/"+dataset+"/FOREST_Forest_"+to_string(DELTA)+".txt";
            ofstream write_index;
            write_index.open(store_forest);
            for(int i=0;i<DELTA;i++)
            {
                write_index<<i<<' '<<forest[i].size()<<endl;
                for(auto &item:forest[i])
                {
                    int node_id=item.ID;
                    FNode node=item;
                    write_index<<node_id<<' '<<node.mu<<' '<<node.father<<' '<<node.childs.size()<<' ';
                    for(auto &child:node.childs)
                    {
                        write_index<<child<<' ';
                    }
                    write_index<<node.cores.size()<<' ';
                    for(auto &v:node.cores)
                    {
                        write_index<<v<<' ';
                    }
                    write_index<<node.non_cores.size()<<' ';
                    for(auto &v:node.non_cores)
                    {
                        write_index<<v<<' ';
                    }
                    write_index<<endl;
                }
            }
            write_index.close();
            cout<<"[STORE INDEX] finish !!!!"<<endl;
        }
        void load_index() override
        {
            cout<<"[LOAD INDEX] start !!!!"<<endl;
            string store_forest="./Index/"+dataset+"/FOREST_Forest_"+to_string(DELTA)+".txt";
            ifstream read_index;
            string line;
            read_index.open(store_forest);
            bool read_node=false;
            int tree_id;
            int tree_len;
            int count=0;
            vertex2node.reserve(DELTA);
            while(getline(read_index, line))
            {
                char* temp;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                if(!read_node)
                {
                    tree_id=atoi(temp);
                    temp = strtok(NULL, " ");
                    tree_len=atoi(temp);
                    if(tree_len) read_node=true;
                    tree_t tree;
                    forest[tree_id]=tree;
                    cout<<"load tree "<<tree_id<<" and it has "<<tree_len<<" nodes "<<endl;
                    count=0;
                    vector<int> v2n;
                    if(tree_len)
                    {
                        for(int i=0;i<n;i++)
                        {
                            v2n.push_back(-1);
                        }
                    }
                    vertex2node.push_back(v2n);
                }
                else
                {
                    FNode temp_node;
                    temp_node.ID=atoi(temp);
                    temp = strtok(NULL, " ");
                    temp_node.mu=atoi(temp);
                    temp = strtok(NULL, " ");
                    temp_node.father=atoi(temp);
                    temp = strtok(NULL, " ");
                    int child_len=atoi(temp);
                    for(int j=0;j<child_len;j++)
                    {
                        temp = strtok(NULL, " ");
                        int child=atoi(temp);
                        temp_node.childs.push_back(child);
                    }
                    temp = strtok(NULL, " ");
                    int core_len=atoi(temp);
                    for(int j=0;j<core_len;j++)
                    {
                        temp = strtok(NULL, " ");
                        int v=atoi(temp);
                        temp_node.cores.push_back(v);
                        vertex2node[tree_id][v]=temp_node.ID;
                    }
                    temp = strtok(NULL, " ");
                    int non_core_len=atoi(temp);
                    for(int j=0;j<non_core_len;j++)
                    {
                        temp = strtok(NULL, " ");
                        int v=atoi(temp);
                        temp_node.non_cores.push_back(v);
                    }
                    forest[tree_id].push_back(temp_node);
                    count++;
                    if(count==tree_len) read_node=false;
                } 
            }
            read_index.close();
            cout<<"[LOAD INDEX] finish!"<<endl;
        }
        void check_query(int epsilon, int mu)
        {
            cout<<"[LOAD INDEX] start loading no and ci!!!!"<<endl;
            string store_no="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            string store_ci="./Index/"+dataset+"/BOTBIN_Cluster_Index.txt";
            ifstream read_index;

            string line;
            read_index.open(store_no);
            while(getline(read_index, line))
            {
                char* temp;
                vector<pair<int,int>> item;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int len=atoi(temp);
                for(int i=0;i<len;i++)
                {
                    temp = strtok(NULL, " ");
                    int sim=atoi(temp);
                    temp = strtok(NULL, " ");
                    int v=atoi(temp);
                    item.push_back(make_pair(sim,v));
                }
                neighbor_order.push_back(item);
            }
            read_index.close();

            vector<vector<pair<int, int>>> bucket;
            read_index.open(store_ci);
            while(getline(read_index, line))
            {
                char* temp;
                vector<pair<int,int>> item;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int len=atoi(temp);
                for(int i=0;i<len;i++)
                {
                    temp = strtok(NULL, " ");
                    int cnt=atoi(temp);
                    temp = strtok(NULL, " ");
                    int v=atoi(temp);
                    item.push_back(make_pair(cnt,v));
                }
                bucket.push_back(item);
            }
            read_index.close();
            cout<<"[LOAD Index] finish loading no and ci !!!!"<<endl;

            cout<<"[CLUSTER] clustering with BOTBIN start !!!"<<endl;
            int locate_bucket=locate(epsilon);
            cout<<"search from bucket["<<locate_bucket<<"]"<<endl;
            vector<int> core_label1;
            for(int i=0;i<n;i++)
            {
                core_label1.push_back(-1);
            }
            unordered_set<int> core_vertice_set;
            for(auto &core:bucket[locate_bucket])
            {
                if(core.first<mu-1) break;
                core_vertice_set.insert(core.second);
            }
            vector<int> watch_vertex={};
            for(auto v:watch_vertex)
            {
                cout<<"vertex "<<v<<" neighbors: ";
                for(auto nv:neighbor_order[v])
                {cout<<"("<<nv.first<<","<<nv.second<<")"<<' ';}
                cout<<endl;
            }
            int cluster_num_of_botbin=0;
            for(auto core:bucket[locate_bucket])
            {
                if(core_label1[core.second]!=-1) continue;
                if(core.first<mu-1) break;
                // cout<<"core: "<<core.second<<endl;
                queue<int> cluster;
                unordered_set<int> cluster_set;
                cluster.push(core.second);
                cluster_set.insert(core.second);
                cluster_num_of_botbin++;
                while(!cluster.empty())
                {
                    int v=cluster.front();
                    cluster.pop();
                    // cout<<v<<" its neighbors: ";
                    for(auto nv:neighbor_order[v])
                    {
                        if(core_label1[nv.second]!=-1) continue;
                        // cout<<"("<<nv.first<<","<<nv.second<<")"<<' ';
                        if(nv.first<epsilon) break;
                        if(core_vertice_set.find(nv.second)!=core_vertice_set.end())//core
                        {
                            if(cluster_set.find(nv.second)==cluster_set.end())
                            {
                                cluster.push(nv.second);
                                cluster_set.insert(nv.second);
                            }
                        }
                        else //non core
                        {

                        }
                    }
                    // cout<<endl;
                }
                int min_id_cluster=n;
                for(auto item:cluster_set)
                {
                    if(item<min_id_cluster) min_id_cluster=item;
                }
                for(auto item:cluster_set)
                {
                    core_label1[item]=min_id_cluster;
                }
                // cout<<endl;
            }
            cout<<"[CLSUTER] finish clustering!!!"<<endl;


            int locate_tree=locate(epsilon);
            cout<<"search from "<<locate_tree<<" th tree"<<endl;
            // printTree(forest[locate_tree]);
            // traverse_print_tree(forest[locate_tree],0);
            vector<int> core_label2(n,-1);
            int cluster_num_of_forest=0;
            traverse_tree_check_query(forest[locate_tree],mu,0,core_label2,cluster_num_of_forest);
            cout<<"botbin find "<<cluster_num_of_botbin<<" clusters"<<endl;
            cout<<"forest find "<<cluster_num_of_forest<<" clusters"<<endl;
            for(int i=0;i<n;i++)
            {
                if(core_label1[i]==-1&&core_label2[i]!=-1)
                {
                    cout<<"vertex "<<i<<" is assign in a cluster of forest but not in a cluster of botbin"<<endl;
                }
                else if(core_label1[i]!=-1&&core_label2[i]==-1)
                {
                    cout<<"vertex "<<i<<" is assign in a cluster of botbin but not in a cluster of forest"<<endl;
                }
                else if(core_label1[i]==-1&&core_label2[i]==-1)
                {
                    continue;
                }
                else if(core_label1[i]!=core_label2[i])
                {
                    cout<<"vertex "<<i<<" belongs to v "<<core_label1[i]<<" of botbin and belongs to v "<<core_label2[i]<<" of forest"<<endl;
                }
                else
                {
                    continue;
                }
            }
            queue<int> cluster;
            unordered_set<int> cluster_set;
            cluster.push(47121);
            cluster_set.insert(47121);
            cout<<"cluster path: "<<endl;
            while(!cluster.empty())
            {
                int v=cluster.front();
                cluster.pop();
                cout<<v<<" its neighbors: ";
                for(auto nv:neighbor_order[v])
                {
                    cout<<"("<<nv.first<<","<<nv.second<<")" ;
                    if(core_label1[nv.second]!=-1) continue;
                    if(nv.first<epsilon) break;
                    cout<<"("<<nv.first<<","<<nv.second<<")" ;
                    if(neighbor_order[nv.second].size()>=mu-1&&neighbor_order[nv.second][mu-2].first>=epsilon) //core
                    {
                        cout<<"is core ";
                        if(cluster_set.find(nv.second)==cluster_set.end())
                        {
                            cluster.push(nv.second);
                            cluster_set.insert(nv.second);
                        }
                    }
                    else //non core
                    {

                    }
                }
                cout<<endl;
            }
        }
        void traverse_tree_check_query(const tree_t&tree, int & mu, int nodeID, vector<int> &core_label2, int & cluster_num)
        {
            if (nodeID >= tree.size()) return;
            const FNode & node=tree[nodeID];
            if(node.mu>=mu)
            {
                unordered_set<int> cluster;
                queue<int> q;
                q.push(node.ID);
                while(!q.empty())
                {
                    int temp_id=q.front();
                    q.pop();
                    const FNode & temp_node=tree[temp_id];
                    for(auto core:temp_node.cores)
                    {cluster.insert(core);} 
                    for(auto child:temp_node.childs)
                    {q.push(child);}
                }
                int min_id_cluster=n;
                for(auto item:cluster)
                {
                    if(item<min_id_cluster) min_id_cluster=item;
                }
                for(auto item:cluster)
                {
                    core_label2[item]=min_id_cluster;
                }
                cluster_num++;
            }
            else
            {
                for(auto child:node.childs)
                {
                    traverse_tree_check_query(tree,mu,child,core_label2,cluster_num);
                }
            }
        }
        void load_similarity_graph()
        {
            string neighbor_order_file="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            ifstream read_index;
            read_index.open(neighbor_order_file);
            string line;
            while(getline(read_index, line))
            {
                char* temp;
                vector<pair<int,int>> item;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int len=atoi(temp);
                for(int i=0;i<len;i++)
                {
                    temp = strtok(NULL, " ");
                    int sim=atoi(temp);
                    temp = strtok(NULL, " ");
                    int v=atoi(temp);
                    item.push_back(make_pair(sim,v));
                }
                neighbor_order.push_back(item);
            }
            read_index.close();
        }

        void print_update_time() override
        {
            if(update_count==0) return;
            cout<<"[PRINT TIME] average update time: "<<total_update_time/double(update_count)<<endl;
        }            
        void collect_clusters_for_update_check(const tree_t &tree, int mu, int node_id, vector<vector<int>> &clusters)
        {
            const FNode &node=tree[node_id];
            if(node.mu>=mu)
            {
                // This node is a cluster frontier at the queried mu. Collect the
                // complete subtree because every descendant belongs to this cluster.
                unordered_set<int> cluster_vertices;
                queue<int> node_queue;
                node_queue.push(node_id);
                while(!node_queue.empty())
                {
                    int current_id=node_queue.front();
                    node_queue.pop();
                    const FNode &current_node=tree[current_id];

                    for(auto vertex:current_node.cores) cluster_vertices.insert(vertex);
                    for(auto vertex:current_node.non_cores) cluster_vertices.insert(vertex);
                    for(auto child:current_node.childs) node_queue.push(child);
                }

                vector<int> cluster(cluster_vertices.begin(),cluster_vertices.end());
                sort(cluster.begin(),cluster.end());
                clusters.push_back(move(cluster));
                return;
            }

            // The current node is below the requested mu; its children may be
            // separate cluster frontiers and must be checked independently.
            for(auto child:node.childs)
            {
                collect_clusters_for_update_check(tree,mu,child,clusters);
            }
        }

        void update_tree(int tree_id, const vector<pair<int,int>> &new_core_v, const vector<int> &new_connectiviy_v, const vector<update_edge_in_graph_t> &update_edges)
        {
            if(new_core_v.empty()&&new_connectiviy_v.empty()&&update_edges.empty()) return;

            int range=tree_id*(PRECISION/DELTA);
            tree_t &tree=forest[tree_id];
            if(tree.empty())
            {
                rebuild_tree(tree_id,new_core_v,new_connectiviy_v,update_edges);
                return;
            }

            if(vertex2node.size()<DELTA) vertex2node.resize(DELTA);
            if(vertex2node[tree_id].size()<n) vertex2node[tree_id].resize(n,-1);

            // Step 1: start a new timestamp round and record the net CT changes.
            // A mark equals forest_update_token only when it belongs to this call,
            // so the O(n) vertex arrays never need to be cleared between updates.
            forest_update_token++;
            if(affected_node_mark.size()<tree.size()) affected_node_mark.resize(tree.size(),0);
            if(affected_root_mark.size()<tree.size()) affected_root_mark.resize(tree.size(),0);
            if(old_root_cache_mark.size()<tree.size()) old_root_cache_mark.resize(tree.size(),0);
            if(old_root_cache.size()<tree.size()) old_root_cache.resize(tree.size(),-1);

            for(auto &changed_vertex:new_core_v)
            {
                int vertex=changed_vertex.first;
                ct_delta_mark[vertex]=forest_update_token;
                update_ct_delta[vertex]=changed_vertex.second;
            }

            vector<int> affected_vertices;
            vector<int> affected_nodes;
            queue<int> vertex_queue;
            queue<int> root_queue;

            // Step 2: initialize the affected-region seeds.
            // CT-changed vertices may move between birth levels. Connectivity
            // vertices remain seeds even when their net CT delta is zero.
            for(auto &changed_vertex:new_core_v)
            {
                int vertex=changed_vertex.first;
                if(affected_vertex_mark[vertex]!=forest_update_token)
                {
                    affected_vertex_mark[vertex]=forest_update_token;
                    affected_vertices.push_back(vertex);
                    vertex_queue.push(vertex);
                }
            }
            for(auto vertex:new_connectiviy_v)
            {
                if(affected_vertex_mark[vertex]!=forest_update_token)
                {
                    affected_vertex_mark[vertex]=forest_update_token;
                    affected_vertices.push_back(vertex);
                    vertex_queue.push(vertex);
                }
            }

            // Step 2.1: also seed endpoints of edges that cross this epsilon.
            // The root-component algorithm does not need to materialize every
            // old/new hierarchy level; its final-graph closure reconstructs all
            // levels inside the selected components in one pass.
            for(auto &edge:update_edges)
            {
                bool old_active=edge.old_similarity>=range;
                bool new_active=edge.new_similarity>=range;
                if(old_active==new_active) continue;

                int edge_vertices[2]={edge.u,edge.v};
                for(int i=0;i<2;i++)
                {
                    int vertex=edge_vertices[i];
                    if(affected_vertex_mark[vertex]!=forest_update_token)
                    {
                        affected_vertex_mark[vertex]=forest_update_token;
                        affected_vertices.push_back(vertex);
                        vertex_queue.push(vertex);
                    }
                }
            }

            // Step 3: compute a closed local reconstruction region.
            // Two alternating expansions are required:
            //   1. old Forest expansion absorbs a complete old root component,
            //      which is necessary because a deletion may split it anywhere;
            //   2. final-graph expansion absorbs components newly connected by an
            //      inserted edge. Newly reached vertices may introduce more old
            //      roots, so the two queues run until both become empty.
            while(!vertex_queue.empty()||!root_queue.empty())
            {
                // Step 3.1: locate the old root of each newly reached vertex.
                // Roots cached while unfolding an old subtree are obtained in O(1);
                // only vertices first discovered through the final graph climb.
                while(!vertex_queue.empty())
                {
                    int vertex=vertex_queue.front();
                    vertex_queue.pop();

                    int node_id=vertex2node[tree_id][vertex];
                    if(node_id>0&&node_id<tree.size())
                    {
                        int root_id=-1;
                        if(old_root_cache_mark[node_id]==forest_update_token)
                        {
                            root_id=old_root_cache[node_id];
                        }
                        else
                        {
                            root_id=node_id;
                            while(tree[root_id].father!=0&&tree[root_id].father!=root_id)
                            {
                                root_id=tree[root_id].father;
                            }
                            old_root_cache_mark[node_id]=forest_update_token;
                            old_root_cache[node_id]=root_id;
                        }

                        if(affected_root_mark[root_id]!=forest_update_token)
                        {
                            affected_root_mark[root_id]=forest_update_token;
                            root_queue.push(root_id);
                        }
                    }

                    // Step 3.2: obtain final CT lazily from old birth mu + delta,
                    // then expand only through final epsilon-edges between cores.
                    int vertex_old_ct=1;
                    node_id=vertex2node[tree_id][vertex];
                    if(node_id>=0&&node_id<tree.size()) vertex_old_ct=tree[node_id].mu;
                    int vertex_new_ct=vertex_old_ct;
                    if(ct_delta_mark[vertex]==forest_update_token)
                    {
                        vertex_new_ct+=update_ct_delta[vertex];
                    }
                    if(vertex_new_ct<2) continue;

                    for(auto &neighbor:neighbor_order[vertex])
                    {
                        if(neighbor.first<range) break;
                        int neighbor_vertex=neighbor.second;

                        int neighbor_old_ct=1;
                        int neighbor_node_id=vertex2node[tree_id][neighbor_vertex];
                        if(neighbor_node_id>=0&&neighbor_node_id<tree.size())
                        {
                            neighbor_old_ct=tree[neighbor_node_id].mu;
                        }
                        int neighbor_new_ct=neighbor_old_ct;
                        if(ct_delta_mark[neighbor_vertex]==forest_update_token)
                        {
                            neighbor_new_ct+=update_ct_delta[neighbor_vertex];
                        }

                        if(neighbor_new_ct<2||affected_vertex_mark[neighbor_vertex]==forest_update_token) continue;
                        affected_vertex_mark[neighbor_vertex]=forest_update_token;
                        affected_vertices.push_back(neighbor_vertex);
                        vertex_queue.push(neighbor_vertex);
                    }
                }

                // Step 3.3: unfold every selected compressed root subtree. Only
                // cores are inserted into the vertex queue; non-core memberships
                // will be regenerated later. Cache root_id for every unfolded node
                // so its cores do not repeat an upward father traversal.
                while(!root_queue.empty())
                {
                    int root_id=root_queue.front();
                    root_queue.pop();
                    queue<int> node_queue;
                    node_queue.push(root_id);
                    while(!node_queue.empty())
                    {
                        int node_id=node_queue.front();
                        node_queue.pop();
                        if(affected_node_mark[node_id]==forest_update_token) continue;
                        affected_node_mark[node_id]=forest_update_token;
                        old_root_cache_mark[node_id]=forest_update_token;
                        old_root_cache[node_id]=root_id;
                        affected_nodes.push_back(node_id);

                        for(auto vertex:tree[node_id].cores)
                        {
                            if(affected_vertex_mark[vertex]!=forest_update_token)
                            {
                                affected_vertex_mark[vertex]=forest_update_token;
                                affected_vertices.push_back(vertex);
                                vertex_queue.push(vertex);
                            }
                        }
                        for(auto child:tree[node_id].childs)
                        {
                            node_queue.push(child);
                        }
                    }
                }
            }

            // Step 4: detach the affected old subtrees from the visible Forest.
            // Unaffected root children remain reachable with their original IDs.
            vector<int> root_childs;
            for(auto child:tree[0].childs)
            {
                if(affected_root_mark[child]!=forest_update_token) root_childs.push_back(child);
            }
            tree[0].childs.swap(root_childs);

            // Step 4.1: clear detached nodes and retain their IDs as a local free
            // list. Reusing these slots prevents tree.size() from growing after
            // every update and keeps all unaffected node IDs unchanged.
            for(auto node_id:affected_nodes)
            {
                tree[node_id].cores.clear();
                tree[node_id].non_cores.clear();
                tree[node_id].childs.clear();
                tree[node_id].father=node_id;
            }

            unordered_map<int,vector<int>> vertices_of_each_ct;
            vector<int> affected_ct_levels;
            vector<int> local_vertices;

            // Step 5: lazily obtain the final CT of affected vertices and place
            // them into final birth levels. CT==1 vertices leave the Forest, but
            // their old vertex2node entries must still be cleared.
            for(auto vertex:affected_vertices)
            {
                int old_ct=1;
                int old_node_id=vertex2node[tree_id][vertex];
                if(old_node_id>=0&&old_node_id<tree.size()) old_ct=tree[old_node_id].mu;
                int new_ct=old_ct;
                if(ct_delta_mark[vertex]==forest_update_token)
                {
                    new_ct+=update_ct_delta[vertex];
                }

                vertex2node[tree_id][vertex]=-1;
                if(new_ct<2) continue;
                if(vertices_of_each_ct.find(new_ct)==vertices_of_each_ct.end())
                {
                    affected_ct_levels.push_back(new_ct);
                }
                vertices_of_each_ct[new_ct].push_back(vertex);
                local_vertices.push_back(vertex);
            }
            sort(affected_ct_levels.begin(),affected_ct_levels.end(),greater<int>());

            // Step 5.1: reset only the local DSU entries. update_uf owns arrays of
            // size n, but Init(local_vertices) changes O(A) entries instead of
            // constructing and initializing a new O(n) UnionFind for every call.
            update_uf.Init(local_vertices);
            vector<int> new_node_ids;
            vector<int> new_root_nodes;
            vector<int> free_node_ids=affected_nodes;

            // Step 6: rebuild the closed region from high mu to low mu.
            // active_vertex_mark identifies cores born at already processed higher
            // levels without clearing an O(n) bitset before each update.
            for(auto mu:affected_ct_levels)
            {
                unordered_map<int,vector<int>> vertex_child_node_map;

                // Step 6.1: before activating this level, locate all adjacent
                // higher-level child components. The current vertex is also a
                // non-core member of each such component at higher query levels.
                for(auto vertex:vertices_of_each_ct[mu])
                {
                    vector<int> child_nodes;
                    for(auto &neighbor:neighbor_order[vertex])
                    {
                        if(neighbor.first<range) break;
                        int neighbor_vertex=neighbor.second;
                        if(active_vertex_mark[neighbor_vertex]!=forest_update_token) continue;

                        int child_node_id=vertex2node[tree_id][neighbor_vertex];
                        tree[child_node_id].non_cores.push_back(vertex);
                        while(tree[child_node_id].father!=0)
                        {
                            child_node_id=tree[child_node_id].father;
                        }
                        child_nodes.push_back(child_node_id);
                    }
                    vertex_child_node_map[vertex]=move(child_nodes);
                }

                // Step 6.2: activate all vertices born at this mu together. This
                // prevents the processing order inside one CT bucket from changing
                // the resulting connected components.
                for(auto vertex:vertices_of_each_ct[mu])
                {
                    active_vertex_mark[vertex]=forest_update_token;
                }

                // Step 6.3: union each newly active vertex with every active
                // epsilon-neighbor, including vertices born in the same bucket.
                for(auto vertex:vertices_of_each_ct[mu])
                {
                    for(auto &neighbor:neighbor_order[vertex])
                    {
                        if(neighbor.first<range) break;
                        int neighbor_vertex=neighbor.second;
                        if(active_vertex_mark[neighbor_vertex]==forest_update_token&&
                           update_uf.Find(vertex)!=update_uf.Find(neighbor_vertex))
                        {
                            update_uf.Unite(vertex,neighbor_vertex);
                        }
                    }
                }

                // Step 6.4: group only the vertices born at this level by their
                // final DSU root. Each group creates one compressed FNode event.
                unordered_map<int,vector<int>> root2vertex;
                for(auto vertex:vertices_of_each_ct[mu])
                {
                    root2vertex[update_uf.Find(vertex)].push_back(vertex);
                }

                // Step 6.5: create the new FNodes. Reuse a detached ID when
                // possible; append only when the new local hierarchy needs more
                // nodes than the old one. Child and vertex IDs remain valid.
                for(auto &component:root2vertex)
                {
                    int node_id;
                    if(!free_node_ids.empty())
                    {
                        node_id=free_node_ids.back();
                        free_node_ids.pop_back();
                    }
                    else
                    {
                        node_id=tree.size();
                    }

                    FNode new_node=FNode(mu,node_id);
                    new_node.father=0;
                    for(auto vertex:component.second)
                    {
                        new_node.cores.push_back(vertex);
                        vertex2node[tree_id][vertex]=node_id;
                        for(auto child:vertex_child_node_map[vertex])
                        {
                            new_node.childs.push_back(child);
                        }
                    }

                    sort(new_node.childs.begin(),new_node.childs.end());
                    new_node.childs.erase(unique(new_node.childs.begin(),new_node.childs.end()),new_node.childs.end());
                    for(auto child:new_node.childs) tree[child].father=node_id;
                    if(node_id<tree.size()) tree[node_id]=move(new_node);
                    else tree.push_back(move(new_node));
                    new_node_ids.push_back(node_id);
                }
            }

            // Step 7: collect only newly rebuilt top components, avoiding a scan
            // over all old and detached FNodes, then splice them under root 0.
            for(auto node_id:new_node_ids)
            {
                if(tree[node_id].father==0)
                {
                    new_root_nodes.push_back(node_id);
                }
            }
            for(auto node_id:new_root_nodes) tree[0].childs.push_back(node_id);

            // Step 8: non-core records may be repeated when several epsilon-edges
            // connect the same lower-level vertex to one higher FNode. Core groups
            // are produced once by DSU and therefore require no sort/unique pass.
            for(auto node_id:new_node_ids)
            {
                sort(tree[node_id].non_cores.begin(),tree[node_id].non_cores.end());
                tree[node_id].non_cores.erase(unique(tree[node_id].non_cores.begin(),tree[node_id].non_cores.end()),tree[node_id].non_cores.end());
            }
        }
        void rebuild_tree(int tree_id, const vector<pair<int,int>> &new_core_v, const vector<int> &new_connectiviy_v, const vector<update_edge_in_graph_t> &update_edges)
        {
            int range=tree_id*(PRECISION/DELTA);
            vector<int> subgraph;
            subgraph.reserve(n);
            vector<bool> vertex_bit_set(n,false);
            vector<int> vertex_node_map(n,-1);
            vector<vector<int>> v_set_of_each_core_threshold(dmax+2);
            UnionFind UF(n);

            // recompute final CT under this epsilon
            for(int vertex=0;vertex<n;vertex++)
            {
                int core_threshold=1;
                for(auto &neighbor:neighbor_order[vertex])
                {
                    if(neighbor.first<range) break;
                    core_threshold++;
                }
                v_set_of_each_core_threshold[core_threshold].push_back(vertex);
            }

            tree_t new_tree;
            FNode root_node=FNode(1,0);
            new_tree.push_back(root_node);
            int node_id=1;

            // rebuild the component hierarchy from high mu to low mu
            for(int mu=dmax+1;mu>1;mu--)
            {
                if(v_set_of_each_core_threshold[mu].empty()) continue;

                unordered_map<int,vector<int>> vertex_child_node_map;
                for(auto vertex:v_set_of_each_core_threshold[mu])
                {
                    vector<int> child_nodes;
                    for(auto &neighbor:neighbor_order[vertex])
                    {
                        if(neighbor.first<range) break;
                        if(!vertex_bit_set[neighbor.second]) continue;

                        int child_node_id=vertex_node_map[neighbor.second];
                        new_tree[child_node_id].non_cores.push_back(vertex);

                        while(new_tree[child_node_id].father!=0)
                        {
                            child_node_id=new_tree[child_node_id].father;
                        }
                        child_nodes.push_back(child_node_id);
                    }
                    vertex_child_node_map[vertex]=child_nodes;
                }

                for(auto vertex:v_set_of_each_core_threshold[mu])
                {
                    subgraph.push_back(vertex);
                    vertex_bit_set[vertex]=true;
                }

                for(auto vertex:v_set_of_each_core_threshold[mu])
                {
                    for(auto &neighbor:neighbor_order[vertex])
                    {
                        if(neighbor.first<range) break;
                        if(vertex_bit_set[neighbor.second]&&
                           UF.Find(vertex)!=UF.Find(neighbor.second))
                        {
                            UF.Unite(vertex,neighbor.second);
                        }
                    }
                }

                unordered_map<int,vector<int>> root2vertex;
                for(auto vertex:v_set_of_each_core_threshold[mu])
                {
                    root2vertex[UF.Find(vertex)].push_back(vertex);
                }

                for(auto &component:root2vertex)
                {
                    FNode new_node=FNode(mu,node_id);
                    new_node.father=0;

                    for(auto vertex:component.second)
                    {
                        new_node.cores.push_back(vertex);
                        vertex_node_map[vertex]=node_id;

                        for(auto child:vertex_child_node_map[vertex])
                        {
                            new_node.childs.push_back(child);
                            new_tree[child].father=node_id;
                        }
                    }

                    new_tree.push_back(new_node);
                    node_id++;
                }
            }

            for(auto &node:new_tree)
            {
                if(node.ID!=0&&node.father==0)
                {
                    new_tree[0].childs.push_back(node.ID);
                }

                sort(node.cores.begin(),node.cores.end());
                node.cores.erase(unique(node.cores.begin(),node.cores.end()),node.cores.end());
                sort(node.non_cores.begin(),node.non_cores.end());
                node.non_cores.erase(unique(node.non_cores.begin(),node.non_cores.end()),node.non_cores.end());
                sort(node.childs.begin(),node.childs.end());
                node.childs.erase(unique(node.childs.begin(),node.childs.end()),node.childs.end());
            }

            forest[tree_id]=move(new_tree);
            if(vertex2node.size()<DELTA)
            {
                vertex2node.resize(DELTA);
            }
            vertex2node[tree_id]=move(vertex_node_map);
        }

        void update_similarity_graph_and_collect_effected_graph(string update_type, int u, int v, vector<tree_update_info_t> &tree_update_info)
        {
            unordered_map<int,vector<int>> vertex_ct_delta;
            vector<update_edge_in_graph_t> similarity_changes;

            unordered_set<int> dirty_neighbor_order;
            dirty_neighbor_order.insert(u);
            dirty_neighbor_order.insert(v);

            // update similarity of (u,w) and (v,w)
            int endpoints[2]={u,v};
            for(int k=0;k<2;k++)
            {
                int endpoint=endpoints[k];
                for(auto &item:neighbor_order[endpoint])
                {
                    int neighbor=item.second;
                    if(neighbor==u||neighbor==v) continue;

                    int old_similarity=item.first;
                    int new_similarity=compute_similarity_exact(endpoint,neighbor);
                    if(old_similarity==new_similarity) continue;

                    item.first=new_similarity;
                    for(auto &reverse_item:neighbor_order[neighbor])
                    {
                        if(reverse_item.second==endpoint)
                        {
                            reverse_item.first=new_similarity;
                            break;
                        }
                    }

                    similarity_changes.push_back({
                        endpoint,
                        neighbor,
                        old_similarity,
                        new_similarity,
                        0
                    });
                    dirty_neighbor_order.insert(neighbor);
                }
            }

            // update similarity edge (u,v)
            if(update_type=="insert")
            {
                int new_uv_similarity=compute_similarity_exact(u,v);
                neighbor_order[u].emplace_back(new_uv_similarity,v);
                neighbor_order[v].emplace_back(new_uv_similarity,u);
                similarity_changes.push_back({
                    u,
                    v,
                    -1,
                    new_uv_similarity,
                    1
                });
            }
            else
            {
                int old_uv_similarity=-1;
                int position_u=-1;
                int position_v=-1;

                for(int i=0;i<neighbor_order[u].size();i++)
                {
                    if(neighbor_order[u][i].second==v)
                    {
                        old_uv_similarity=neighbor_order[u][i].first;
                        position_u=i;
                        break;
                    }
                }

                for(int i=0;i<neighbor_order[v].size();i++)
                {
                    if(neighbor_order[v][i].second==u)
                    {
                        position_v=i;
                        break;
                    }
                }

                neighbor_order[u].erase(neighbor_order[u].begin()+position_u);
                neighbor_order[v].erase(neighbor_order[v].begin()+position_v);

                similarity_changes.push_back({
                    u,
                    v,
                    old_uv_similarity,
                    -1,
                    -1
                });
            }

            // restore neighbor order
            for(auto vertex:dirty_neighbor_order)
            {
                sort(neighbor_order[vertex].begin(),neighbor_order[vertex].end(),cmp_order_forest);
            }

            // Collect only changes that cross at least one discretized epsilon.
            // Similarities staying in the same interval have already been updated
            // in neighbor_order, but they cannot change any Forest tree.
            int bucket_width=PRECISION/DELTA;
            for(auto &edge:similarity_changes)
            {
                int old_max_tree=0;
                int new_max_tree=0;

                if(edge.old_similarity>=0)
                {
                    old_max_tree=edge.old_similarity/bucket_width;
                    old_max_tree=min(old_max_tree,DELTA-1);
                }
                if(edge.new_similarity>=0)
                {
                    new_max_tree=edge.new_similarity/bucket_width;
                    new_max_tree=min(new_max_tree,DELTA-1);
                }

                if(old_max_tree==new_max_tree) continue;

                int change=1;
                int begin_tree=old_max_tree+1;
                int end_tree=new_max_tree;
                if(old_max_tree>new_max_tree)
                {
                    change=-1;
                    begin_tree=new_max_tree+1;
                    end_tree=old_max_tree;
                }

                // The same similarity edge may affect several consecutive trees.
                // Store one copy in each tree so later reconstruction receives
                // exactly the edges whose active state changed at its epsilon.
                for(int tree_id=begin_tree;tree_id<=end_tree;tree_id++)
                {
                    tree_update_info[tree_id].update_edges.push_back(edge);
                }

                int edge_vertices[2]={edge.u,edge.v};
                for(int i=0;i<2;i++)
                {
                    int vertex=edge_vertices[i];
                    if(vertex_ct_delta.find(vertex)==vertex_ct_delta.end())
                    {
                        vertex_ct_delta[vertex]=vector<int>(DELTA,INT_MAX);
                    }

                    for(int tree_id=begin_tree;tree_id<=end_tree;tree_id++)
                    {
                        if(vertex_ct_delta[vertex][tree_id]==INT_MAX)
                        {
                            vertex_ct_delta[vertex][tree_id]=change;
                        }
                        else
                        {
                            vertex_ct_delta[vertex][tree_id]+=change;
                        }
                    }
                }
            }

            // Convert the per-vertex delta arrays into the three tree inputs.
            // INT_MAX means untouched; zero means connectivity changed although
            // insertions and deletions cancelled in the final CT value.
            for(auto &vertex_delta:vertex_ct_delta)
            {
                int vertex=vertex_delta.first;
                for(int tree_id=1;tree_id<DELTA;tree_id++)
                {
                    int ct_delta=vertex_delta.second[tree_id];
                    if(ct_delta==INT_MAX) continue;

                    tree_update_info[tree_id].new_connectiviy_v.push_back(vertex);
                    if(ct_delta!=0)
                    {
                        tree_update_info[tree_id].new_core_v.push_back(make_pair(vertex,ct_delta));
                    }
                }
            }
        }
        
        void update(string update_type, int u, int v, int update_way) override
        {
            // Load the old persistent state only before the first update. Later
            // calls continue from the neighbor_order and Forest already updated
            // in memory by the preceding edge operations.
            
            if(neighbor_order.empty())
            {
                cout<<"[UPDATE] load similarity graph"<<endl;
                load_similarity_graph();
            }
            if(update_way!=0&&vertex2node.empty())
            {
                cout<<"[UPDATE] load similarity graph"<<endl;
                load_index();
            }
            cout<<"[UPDATE] update similairty graph"<<endl;
            vector<tree_update_info_t> tree_update_info(DELTA);
            update_similarity_graph_and_collect_effected_graph(
                update_type,
                u,
                v,
                tree_update_info);
            
            cout<<"[UPDATE] now update index begin!!!"<<endl;
            auto update_begin=std::chrono::steady_clock().now();
            if(update_way==0)
            {//rebuild whole index
                // Reconstruct the complete Forest from final neighbor_order.
                // construct_forest also refreshes every tree's vertex2node map,
                // so later edge updates can continue from this rebuilt state.
                construct_forest();
            }
            else if(update_way==1)
            {//rebuild affected tree 
                for(int tree_id=1;tree_id<DELTA;tree_id++)
                {
                    tree_update_info_t &update_info=tree_update_info[tree_id];
                    if(update_info.new_connectiviy_v.empty()) continue;

                    rebuild_tree(
                        tree_id,
                        update_info.new_core_v,
                        update_info.new_connectiviy_v,
                        update_info.update_edges);
                }
            }
#if 0
            // Component-level maintenance and its verification mode are kept
            // for future work. Current experiments expose only way 0 and 1.
            else if(update_way==2)
            {//update affected tree
                for(int tree_id=1;tree_id<DELTA;tree_id++)
                {
                    tree_update_info_t &update_info=tree_update_info[tree_id];
                    if(update_info.new_connectiviy_v.empty()) continue;

                    update_tree(
                        tree_id,
                        update_info.new_core_v,
                        update_info.new_connectiviy_v,
                        update_info.update_edges);
                }
            }
            else if(update_way==3)
            {//update affected tree and verify it by full reconstruction
                for(int tree_id=1;tree_id<DELTA;tree_id++)
                {
                    tree_update_info_t &update_info=tree_update_info[tree_id];
                    if(update_info.new_connectiviy_v.empty()) continue;

                    // Step 1: produce the incremental result that should remain as
                    // the in-memory state for the next edge update.
                    update_tree(
                        tree_id,
                        update_info.new_core_v,
                        update_info.new_connectiviy_v,
                        update_info.update_edges);

                    tree_t incremental_tree=forest[tree_id];
                    vector<int> incremental_vertex2node=vertex2node[tree_id];

                    // Step 2: independently rebuild the same tree from the final
                    // neighbor_order. This result is used only as the oracle.
                    rebuild_tree(
                        tree_id,
                        update_info.new_core_v,
                        update_info.new_connectiviy_v,
                        update_info.update_edges);

                    bool correct=true;
                    int wrong_mu=-1;
                    vector<vector<int>> incremental_clusters;
                    vector<vector<int>> rebuilt_clusters;

                    // Step 3: compare semantic clusters for every valid mu. Node
                    // IDs and child order are deliberately ignored because two
                    // equivalent compressed trees may use different storage IDs.
                    for(int mu=2;mu<=dmax+1;mu++)
                    {
                        incremental_clusters.clear();
                        rebuilt_clusters.clear();

                        collect_clusters_for_update_check(
                            incremental_tree,
                            mu,
                            0,
                            incremental_clusters);
                        collect_clusters_for_update_check(
                            forest[tree_id],
                            mu,
                            0,
                            rebuilt_clusters);

                        sort(incremental_clusters.begin(),incremental_clusters.end());
                        sort(rebuilt_clusters.begin(),rebuilt_clusters.end());
                        if(incremental_clusters!=rebuilt_clusters)
                        {
                            correct=false;
                            wrong_mu=mu;
                            break;
                        }
                    }

                    if(correct)
                    {
                        cout<<"[UPDATE CHECK] tree "<<tree_id<<" passed"<<endl;
                    }
                    else
                    {
                        cout<<"[UPDATE CHECK ERROR] tree "<<tree_id
                            <<" mu "<<wrong_mu
                            <<" incremental clusters "<<incremental_clusters.size()
                            <<" rebuilt clusters "<<rebuilt_clusters.size()<<endl;
                    }

                    // Step 4: restore the incremental result after checking. The
                    // next update therefore verifies a genuinely continuous chain
                    // of local updates instead of continuing from rebuilt trees.
                    forest[tree_id]=move(incremental_tree);
                    vertex2node[tree_id]=move(incremental_vertex2node);
                }
            }
#endif
            else
            {
                cout<<"[ERROR] no such update way !!!"<<endl;
                return;
            }
            auto update_end=std::chrono::steady_clock().now();
            cout<<"[UPDATE] update index finish !!!"<<endl;
            tms update_time=update_end-update_begin;
            cout<<"[UPDATE] update time : "<<update_time.count()<<endl;
            total_update_time+=update_time.count();
            update_count++;
        }

};
