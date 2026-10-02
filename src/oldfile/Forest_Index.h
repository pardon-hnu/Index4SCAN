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
#include "global.h"
#include "Index.h"
#include "RandList.h"
#include "SDSU.h"
using namespace std;
// #define DEBUG 
// #define PRINT 
#define DEBUG_CLUSTERING
#define NEW_BUILD
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
        void construct_tree(int tree_id)
        {
            //range: [1-i/delta, 1-(i+1)/delta]
            double range=1.00-double(tree_id+1)/double(DELTA);

            //init i-th tree
            tree_t tree;
            FNode root_node=FNode(1,0);
            tree.push_back(root_node);
            int ID=1;

            //init a map: vertex id->tree node id
            // unordered_map<int,int> vertex_node_map;
            vector<int> vertex_node_map(n,-1);

            //init a set: node without father
            unordered_set<int> rest_node;

            //get range-neighbors number for all vertices
            unordered_map<int,vector<int>> Vset=dsu.mu_decomposition_and_update(neighbor_order,range);

            //bottom up
            for(int temp_mu=dmax+1;temp_mu>1;--temp_mu)
            {
                //ensure exists (epsilon,i)-level
                if(Vset.find(temp_mu)==Vset.end()) continue;
                #ifdef PRINT
                cout<<"[PROCESS] now deal "<<temp_mu<<" level "<<endl;
                #endif
                vector<int> V_mu=Vset[temp_mu];
                #ifdef DEBUG
                cout<<"V mu: ";
                for(auto &v:V_mu)
                {
                    cout<<v<<' ';
                }
                cout<<endl;
                #endif
                /*
                    get the child node for each vertex v in V_mu
                    if v' neighbor nv is in next-level tree node
                    find its corresponding node's root vertex
                */
                unordered_map<int,vector<int>> vertex_subnode_map;
                for(auto &v:V_mu)
                {
                    vector<int> subnode;
                    for(auto &nv:neighbor_order[v])
                    {
                        //if nv is in the next-level 
                        if(dsu.nodes[nv.second].cur_mu>dsu.nodes[v].cur_mu&&nv.first>=range)
                        {
                            int ru=dsu.find(dsu.nodes[nv.second]);
                            int hook_ru=dsu.nodes[ru].hook;
                            int node_ru=vertex_node_map[hook_ru];
                            subnode.push_back(node_ru);
                        }
                    }
                    vertex_subnode_map[v]=subnode;
                }
                #ifdef DEBUG
                cout<<"[VERTEX_TO_SUBNODE] print "<<endl;
                for(auto &v:V_mu)
                {
                    cout<<"v: "<<v<<" subnodes: ";
                    for(auto &subnode_id:vertex_subnode_map[v])
                    {
                        cout<<subnode_id<<' ';
                    }
                    cout<<endl;
                }
                #endif
                /*
                    if a componet is valid in (epsilon+0.1,i)-clusters
                    it must be valid in (epsilon,i)-clusters
                    such connectivity checking is redundant
                */
                vector<int> re_V;
                vector<int> ne_V;
                if(tree_id==0)
                {
                    /*
                        No redundant computation when constructing first tree
                        then we need to init the union-find for all V_mu
                    */
                    for(auto &v:V_mu)
                    {
                        dsu.make_set(v);
                    }
                    ne_V=V_mu;
                }
                else
                {
                    /*
                        checked connectivity in (>epsilon,>mu)-clusters 
                        need not to be checked again
                    */
                    for(auto &v:V_mu)
                    {
                        if(dsu.nodes[v].pre_mu>=temp_mu)
                        {
                            dsu.nodes[v].save_componet_and_init();
                            re_V.push_back(v);
                        }
                        else
                        {
                            dsu.nodes[v].init();
                            ne_V.push_back(v);
                        }
                    }
                }
                #ifdef DEBUG
                cout<<"redundant V: ";
                for(auto &v:re_V)
                    cout<<v<<' ';
                cout<<endl;
                for(int temp_ite=dmax+1;temp_ite>=temp_mu;--temp_ite)
                {
                    if(Vset.find(temp_ite)==Vset.end()) continue;
                    for(auto &v:Vset[temp_ite])
                    {
                        cout<<"vertex id: "<<v<<" root: "<<dsu.nodes[v].father<<" rank: "<<
                        dsu.nodes[v].rank<<" componet: "<<dsu.nodes[v].componet<< 
                            " hook: "<<dsu.nodes[v].hook<<" cur_mu: "<<dsu.nodes[v].cur_mu
                            <<" pre_mu: "<<dsu.nodes[v].pre_mu<<endl;
                    }
                }
                #endif
                /*
                    connect vertex v and its subnode's core vertex
                */
                // for(auto &v:V_mu)
                // {
                //     for(auto &subnode_ID:vertex_subnode_map[v])
                //     {
                //         for(auto &sv:tree[subnode_ID].cores)
                //         {
                //             dsu.unite(v,sv);
                //             #ifdef DEBUG
                //             cout<<"unite vertex "<<v<<" and vertex "<<sv<<endl;
                //             #endif
                //         }
                //     }
                // }
                /*
                    deal with each vertex v in the new V
                    if its neighbor u with higher mu 
                    unite u and v, the are in the same tree (maybe not same level)
                */
                for(auto &v:ne_V)
                {
                    for(auto &nv:neighbor_order[v])
                    {
                        if(dsu.nodes[nv.second].cur_mu>=dsu.nodes[v].cur_mu&&nv.first>=range)
                        {
                            dsu.unite(v,nv.second);
                            #ifdef DEBUG
                            cout<<"union "<<v<<" and "<<nv.second<<endl;
                            #endif
                        }
                    }
                }
                /*
                    deal with each vertex v in the redundant V
                    unite v and its componet which is united in (>epsilon,>mu)
                */
                for(auto &v:re_V)
                {
                    dsu.unite(v,dsu.nodes[v].componet);
                    #ifdef DEBUG
                    cout<<"union "<<v<<" and "<<dsu.nodes[v].componet<<endl;
                    #endif
                }
                /*
                    creat tree node, link, and hook 
                    (1) each node corresponding to a root id
                    (2) add v in V mu into node's v_set
                    (3) set child and father of this node according to vertex_subnode_map
                */
                #ifdef DEBUG
                cout<<"After UNION!!!!!!!"<<endl;
                for(int temp_ite=dmax+1;temp_ite>=temp_mu;--temp_ite)
                {
                    if(Vset.find(temp_ite)==Vset.end()) continue;
                    for(auto &v:Vset[temp_ite])
                    {
                        cout<<"vertex id: "<<v<<" root: "<<dsu.nodes[v].father<<" rank: "<<
                        dsu.nodes[v].rank<<" componet: "<<dsu.nodes[v].componet<< 
                            " hook: "<<dsu.nodes[v].hook<<" cur_mu: "<<dsu.nodes[v].cur_mu
                            <<" pre_mu: "<<dsu.nodes[v].pre_mu<<endl;
                    }
                }
                #endif

                for(auto &v:V_mu)
                {
                    int rv=dsu.find(v);
                }

                unordered_map<int,int> temp_root_node_map;
                for(auto &v:V_mu)
                {
                    int rv=dsu.find(dsu.nodes[v]);
                    if(temp_root_node_map.find(rv)==temp_root_node_map.end())
                    {
                        FNode temp_node=FNode(temp_mu,ID);
                        temp_root_node_map[rv]=ID;
                        rest_node.insert(ID);
                        tree.push_back(temp_node);
                        #ifdef DEBUG
                        cout<<"creat a node with ID "<<ID<<endl;
                        #endif
                        ID++;
                    }
                    int rv_node_id=temp_root_node_map[rv];
                    tree[rv_node_id].cores.push_back(v);
                    vertex_node_map[v]=rv_node_id;
                    #ifdef DEBUG
                    cout<<"add vertex "<<v<<" into v set of node "<<rv_node_id<<endl;
                    #endif
                    for(auto &subnode_ID:vertex_subnode_map[v])
                    {
                        if(tree[subnode_ID].father!=subnode_ID) continue;
                        tree[rv_node_id].childs.push_back(subnode_ID);
                        tree[subnode_ID].father=rv_node_id;
                        
                        #ifdef DEBUG
                        cout<<"add node "<<subnode_ID<<" into child set of node "<<rv_node_id<<endl;
                        cout<<"set node "<<rv_node_id<<" as father of node "<<subnode_ID<<endl;
                        #endif
                        if(rest_node.find(subnode_ID)!=rest_node.end()) rest_node.erase(subnode_ID);
                    }
                }
                
                #ifdef DEBUG
                cout<<"state after create node"<<endl;
                for(int temp_ite=dmax+1;temp_ite>=temp_mu;--temp_ite)
                {
                    if(Vset.find(temp_ite)==Vset.end()) continue;
                    for(auto &v:Vset[temp_ite])
                    {
                        cout<<"vertex id: "<<v<<" root: "<<dsu.nodes[v].father<<" rank: "<<
                        dsu.nodes[v].rank<<" componet: "<<dsu.nodes[v].componet<< 
                            " hook: "<<dsu.nodes[v].hook<<" cur_mu: "<<dsu.nodes[v].cur_mu
                            <<" pre_mu: "<<dsu.nodes[v].pre_mu<<endl;
                    }
                }
                #endif

                //getting non core for each node in this level
                // for(auto &item:temp_root_node_map)
                // {
                    // unordered_set<int> sub_non_cores;
                    // int rootID=item.second;
                    // for(auto &childID:tree[rootID].childs)
                    // {
                    //     for(auto &nc:tree[childID].non_cores)
                    //     {
                    //         sub_non_cores.insert(nc);
                    //     }
                    // }
                    // for(auto &v:tree[rootID].cores)
                    // {
                    //     for(auto &nv:neighbor_order[v])
                    //     {
                    //         if(nv.first>=range&&dsu.nodes[nv.second].cur_mu<dsu.nodes[v].cur_mu&&sub_non_cores.find(nv.second)==sub_non_cores.end())
                    //         {
                    //             tree[rootID].non_cores.insert(nv.second);
                    //         }
                    //     }
                    // }
                // }
                //update DSU
                dsu.update_dsu(V_mu);

                #ifdef DEBUG
                cout<<"state of vertex in higher node"<<endl;
                for(int temp_ite=dmax+1;temp_ite>=temp_mu;--temp_ite)
                {
                    if(Vset.find(temp_ite)==Vset.end()) continue;
                    for(auto &v:Vset[temp_ite])
                    {
                        cout<<"vertex id: "<<v<<" root: "<<dsu.nodes[v].father<<" rank: "<<
                        dsu.nodes[v].rank<<" componet: "<<dsu.nodes[v].componet<< 
                            " hook: "<<dsu.nodes[v].hook<<" cur_mu: "<<dsu.nodes[v].cur_mu
                            <<" pre_mu: "<<dsu.nodes[v].pre_mu<<endl;
                    }
                }
                #endif
            }
            for(int i=0;i<n;i++)
            {
                // int max_mu_id=-1;
                // int max_mu=-1;
                for(int j=0;j<neighbor_order[i].size();j++)
                {     
                    if(neighbor_order[i][j].first<range) break;
                    int nei_id=neighbor_order[i][j].second;
                    if(dsu.nodes[i].cur_mu<dsu.nodes[nei_id].cur_mu)
                    {
                        int relate_node_id=vertex_node_map[nei_id];
                        tree[relate_node_id].non_cores.push_back(i);   
                    }
                    // if(max_mu<dsu.nodes[nei_id].cur_mu&&dsu.nodes[i].cur_mu<dsu.nodes[nei_id].cur_mu)
                    // {
                    //     max_mu=dsu.nodes[nei_id].cur_mu;
                    //     max_mu_id=nei_id;
                    // }
                }
                // if(max_mu>1)
                // {
                //     int relate_node_id=vertex_node_map[max_mu_id];
                //     tree[relate_node_id].non_cores.insert(i);
                // }
            }
            /*
                create root node of this tree
                link nodes without father to root node
            */
            #ifdef PRINT
            cout<<"[PROCESS] now deal the rest node"<<endl;
            #endif
            
            for(auto &rest_node_id:rest_node)
            {
                tree[rest_node_id].father=0;
                tree[0].childs.push_back(rest_node_id);
                #ifdef DEBUG
                cout<<"add node "<<rest_node_id<<" into child set of node 0 "<<endl;
                cout<<"set node 0 as father of node "<<rest_node_id<<endl;
                #endif
            }
            for(auto item:tree)
            {
                sort(item.cores.begin(), item.cores.end());
                item.cores.erase(unique(item.cores.begin(), item.cores.end()), item.cores.end());
                sort(item.non_cores.begin(), item.non_cores.end());
                item.non_cores.erase(unique(item.non_cores.begin(), item.non_cores.end()), item.non_cores.end());
                sort(item.childs.begin(), item.childs.end());
                item.childs.erase(unique(item.childs.begin(), item.childs.end()), item.childs.end());
            }

            forest[tree_id]=tree;
        }  
        void construct_forest2()
        {
            cout<<"[CONSTRUCT FOREST] start"<<endl;
            auto construct_forest_begin=std::chrono::steady_clock().now();
            vector<int> subgraph;
            subgraph.reserve(n);
            vector<bool> vertex_bit_set;
            vertex_bit_set.resize(n);
            std::fill_n(vertex_bit_set.begin(), vertex_bit_set.size(), false);
            UnionFind UF(n);
            for(int tree_id=0;tree_id<DELTA-1;tree_id++)
            {
                cout<<"[CONSTRUCT FOREST] now construct "<<tree_id<<" th tree"<<endl;
                 //range: [1-i/delta, 1-(i+1)/delta]
                double range=1.00-double(tree_id+1)/double(DELTA);

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
        void construct_forest()
        {
            cout<<"[CONSTRUCT FOREST] start"<<endl;
            auto construct_forest_begin=std::chrono::steady_clock().now();
            for(int i=0;i<DELTA-1;i++)
            {
                cout<<"[CONSTRUCT FOREST] now construct "<<i<<" th tree"<<endl;
                construct_tree(i);
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
            #ifdef NEW_BUILD
                construct_forest2();
            #else
                construct_forest();
            #endif
        }
        int locate(float epsilon)
        {
            if(epsilon==0.9f) return 0;
            double num=DELTA;
            int i=1;
            for(i=DELTA;i>0;i--)
            {
                double similairty=1.0-double(i)/num;
                if(epsilon<similairty)
                {
                    break;
                }
            }
            return i;
        }
        void traverse_tree_query(const tree_t&tree, int & mu, bool & clusterFLAG, int nodeID, int & clusterID)
        {
            if (nodeID >= tree.size()) return;
            const FNode & node=tree[nodeID];
            bool init_flag=false;
            #ifdef DEBUG
                cout<<"node id: "<<node.ID<<" mu: "<<node.mu<<" cluster flag: "<<clusterFLAG<<" cluster id"<<clusterID<<endl;
            #endif
            traverse_node_num++;
            if(node.mu>=mu)
            {
                if(!clusterFLAG)
                {
                    clusterID++;
                    clusterFLAG=true;
                    init_flag=true;
                    #ifdef DEBUG
                        cout<<"--begin clustering, cluster id: "<<clusterID<<endl;
                    #endif
                }
                for(auto &v:node.cores)
                {
                    traverse_core_num++;
                    labels[v]=clusterID;
                }
                for(auto &v:node.non_cores)
                {
                    traverse_non_core_num++;
                    labels[v]=clusterID;
                }
            }
            for(auto &childID:node.childs)
            {
                traverse_tree_query(tree,mu,clusterFLAG,childID,clusterID);
            }
            if(init_flag)
            {
                clusterFLAG=false;
                #ifdef DEBUG
                    cout<<"--finish clustering, cluster id: "<<clusterID<<endl;
                    for(int i=0;i<n;i++)
                    {
                        cout<<"(node: "<<i<<",label: "<<labels[i]<<") ";
                    }
                    cout<<endl;
                #endif
            }
        }
        void forest_query(float epsilon, int mu)
        {
            clear_cluster_result();
            traverse_node_num=0;
            traverse_core_num=0;
            traverse_non_core_num=0;
            cout<<"[QUERY] Start"<<endl;
            auto query_begin=std::chrono::steady_clock().now();
            int locate_tree=locate(epsilon);
            cout<<"search from "<<locate_tree<<" th tree"<<endl;
            bool flag1=false;
            bool flag2=false;
            int clusterID=-1;
            traverse_tree_query(forest[locate_tree],mu,flag1,0,clusterID);
            auto query_end=std::chrono::steady_clock().now();
            tms query_time=query_end-query_begin;
            cout<<"[QUERY] finish"<<endl;
            cout<<"[TIME COST] "<<query_time.count()<<" s."<<endl;
            total_cluster_time+=query_time.count();
            query_count++;
            #ifdef EVALUATION_CLUSTERING_QUALITY
                evaluate_clustering_quality(epsilon,mu);
            #endif
            cout<<"[INFORMATION DEBUG] traverse node num: "<<traverse_node_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse core num: "<<traverse_core_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse non core num: "<<traverse_non_core_num<<endl;
        }
        void query(float epsilon, int mu) override
        {
            #ifdef DEBUG_CLUSTERING
                check_query(epsilon,mu);
            #else
                forest_query(epsilon,mu);
            #endif
        }
        void evaluate_clustering_quality(float epsilon, int mu) 
        {
            string groundtruth_file="Result/labels/"+dataset+"/GS-Index-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".txt";
            vector<int> groundtruth;
            groundtruth.reserve(n);
            ifstream read_groundtruth;
            read_groundtruth.open(groundtruth_file);
            string line;
            getline(read_groundtruth, line);
            char* temp;
            temp = strtok(const_cast<char*>(line.c_str()), " ");
            while (temp != NULL) {
                int label = atoi(temp);
                groundtruth.push_back(label);
                temp = strtok(NULL, " ");
            }
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
                cout<<"print the "<<i<<" th tree corresponding range ("<<1.00-float(i)/float(DELTA)<<","<<1.00-float(i+1)/float(DELTA)<<"]"<<endl;
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
            // string store_rnk="./Index/"+dataset+"/BOTBIN_GRAPH_HASH.txt";
            // string store_no="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            string store_forest="./Index/"+dataset+"/FOREST_Forest_"+to_string(DELTA)+".txt";
            ifstream read_index;

            // read_index.open(store_rnk);
            string line;
            // while(getline(read_index, line))
            // {
            //     char* temp;
            //     vector<pair<int,float>> item;
            //     temp = strtok(const_cast<char*>(line.c_str()), " ");
            //     int hash_v=atoi(temp);
            //     rnk.push_back(hash_v);
            // }
            // int* rnk_id = new int[n];
            // for (int i = 0; i < n; i++) {
            //     rnk_id[i] = i;
            // }
            // sort(rnk_id, rnk_id + n, cmp_rnk);
            // sketch.resize(n);

            // for (int i = 0; i < n; i++) {

            //     int u = rnk_id[i];
            //     if (sketch[u].size() < k)
            //         sketch[u].push_back(i);

            //     for (int j = 0; j < degree[u]; j++) {
            //         int v = graph[u][j];
            //         if (sketch[v].size() >= k) {
            //             continue;
            //         }
            //         sketch[v].push_back(i);
            //     }
            // }
            // for(int i=0;i<n;i++)
            // {
            //     sort(sketch[i].begin(),sketch[i].end());
            // }
            // read_index.close();

            // read_index.open(store_no);
            // while(getline(read_index, line))
            // {
            //     char* temp;
            //     vector<pair<float,int>> item;
            //     temp = strtok(const_cast<char*>(line.c_str()), " ");
            //     int len=atoi(temp);
            //     for(int i=0;i<len;i++)
            //     {
            //         temp = strtok(NULL, " ");
            //         float sim=atof(temp);
            //         temp = strtok(NULL, " ");
            //         int v=atoi(temp);
            //         item.push_back(make_pair(sim,v));
            //     }
            //     neighbor_order.push_back(item);
            // }
            // read_index.close();

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
                    read_node=true;
                    count=0;
                    // cout<<"read tree: "<<tree_id<<" and number of nodes: "<<tree_len<<endl;
                }
                else
                {
                    if(count==0)
                    {
                        tree_t tree;
                        forest[tree_id]=tree;
                    }
                    if(count<tree_len)
                    {
                        count++;
                    }
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
                    // cout<<"read node: "<<temp_node.ID<<" mu: "<<temp_node.mu<<" father: "<<temp_node.father<<" childs: ";
                    // for(auto &c:temp_node.childs)
                    //     cout<<c<<" ";
                    // cout<<"cores: ";
                    // for(auto &v:temp_node.cores)
                    //     cout<<v<<" ";
                    // cout<<endl;
                    forest[tree_id].push_back(temp_node);
                    if(count==tree_len) read_node=false;
                } 
            }
            read_index.close();
            // cout<<"--------------------------------------------"<<endl;
            // print_index();
            cout<<"[LOAD INDEX] finish!"<<endl;
        }
        void check_query(float epsilon, int mu)
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
                vector<pair<float,int>> item;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int len=atoi(temp);
                for(int i=0;i<len;i++)
                {
                    temp = strtok(NULL, " ");
                    float sim=atof(temp);
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
            vector<int> watch_vertex={207419, 307677, 114700};
            for(auto v:watch_vertex)
            {
                cout<<"vertex "<<v<<" neighbors: ";
                for(auto nv:neighbor_order[v])
                {cout<<"("<<nv.first<<","<<nv.second<<")"<<' ';}
                cout<<endl;
            }
            for(auto core:bucket[locate_bucket])
            {
                if(core_label1[core.second]!=-1) continue;
                if(core.first<mu-1) break;
                // cout<<"core: "<<core.second<<endl;
                queue<int> cluster;
                unordered_set<int> cluster_set;
                cluster.push(core.second);
                cluster_set.insert(core.second);
                int temp_cluster_number=0;
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
            printTree(forest[locate_tree]);
            traverse_print_tree(forest[locate_tree],0);
            vector<int> core_label2(n,-1);
            traverse_tree_check_query(forest[locate_tree],mu,0,core_label2);
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
        }
        void traverse_tree_check_query(const tree_t&tree, int & mu, int nodeID, vector<int> &core_label2)
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
            }
            else
            {
                for(auto child:node.childs)
                {
                    traverse_tree_check_query(tree,mu,child,core_label2);
                }
            }
        }
};