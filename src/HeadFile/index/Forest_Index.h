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

        vector<vector<pair<int, int>>> neighbor_order;
        vector<tree_t> forest;

        int total_node_num=0;
        int traverse_node_num=0;
        int traverse_core_num=0;
        int traverse_non_core_num=0;

        double core_vertices_space=0.0;
        double non_core_vertices_space=0.0;

        double total_ari=0.0;
        int query_count=0;
    public:
        SCAN_Forest_Index(string dataset,int DELTA=10)
        {
            this->DELTA=DELTA;
            forest.resize(DELTA);
            init_graph(dataset);
            k=(int)(ceil((1/(2*pow(RHO,2)))*(log(2/FAILURE_PB))));
            dsu=SCAN_DSU(n);
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
                forest[tree_id]=tree;
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
};