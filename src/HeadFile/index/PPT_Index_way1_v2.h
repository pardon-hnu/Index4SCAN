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
#include <index/PPT.h>
#include <util/RandList.h>
#include <util/RStarTree.h>
// #define CHECK_RESULT
// #define TEST_STRUCTURE
// #define TEST_CORRECTNESS
typedef RStarTree<int, 2, 32, 64> 			RTree;
typedef RTree::BoundingBox BoundingBox;
BoundingBox bounds(int x, int y, int w, int h)
{
	BoundingBox bb;
	
	bb.edges[0].first  = x;
	bb.edges[0].second = x + w;
	
	bb.edges[1].first  = y;
	bb.edges[1].second = y + h;
	
	return bb;
}
struct Visitor {
	int count;
	bool ContinueVisiting;
	
	Visitor() : count(0), ContinueVisiting(true) {};
	
	void operator()(const RTree::Leaf * const leaf) 
	{
		count++;
	}
};

struct ppt_update_edge_t
{
    int u;
    int v;
    int old_similarity;
    int new_similarity;
};

bool cmp_order_space_index(pair<float,int> x,pair<float,int> y)
{
	return x.first>y.first;
}
class SCAN_PPT_Index : public SCAN_Index
{
    private:
        int k=0;
        float RHO=0.1;
        float FAILURE_PB=0.001;
        int DELTA=10;
        string format;

        static vector<int> rnk;
        vector<vector<int> > sketch;

        vector<vector<pair<int, int>>> neighbor_order;
        vector<vector<int>> interval_num;

        vector<PPT> ppts;
        vector<vector<int>> cluster_slices;
        vector<vector<int>> non_cores;
        int ppt_max_id=0;
        unordered_map<pair<int,int>,int,PairHash> point2ppt;
        vector<vector<int>> vertex2ppt;
        vector<vector<int>> vertex2cs;
        vector<int> cs2ppt;
        vector<int> free_update_cs_ids;
        
        int * max_mu_under_epsilon=nullptr; // when fix a certain epsilon, the max mu that (epsilon,mu) is in space.
        vector<vector<int>> ppt_matrix;  // PPT[i][j]!=-1 :  there exist a PPT(i,j+2)

        vector<vector<int>> table_ppts;
        RTree r_star_tree;

        int traverse_node_num=0;
        int traverse_core_num=0;
        int traverse_non_core_num=0;
        int traverse_edge_num=0;
        int traverse_unite_num=0;
        int traverse_pc_find_num=0;
        float traverse_unite_time=0.0;
        float traverse_pc_find_time=0.0;
        double core_vertices_space=0.0;
        double non_core_vertices_space=0.0;
        double edge_space=0.0;
        double connectivity_space=0.0;
        double total_ari=0.0;
        double total_collection_ppt_time=0.0;
        int query_count=0;
        int update_count=0;
        double update_profile_select_time=0.0;
        double update_profile_cleanup_time=0.0;
        double update_profile_state_time=0.0;
        double update_profile_uf1_time=0.0;
        double update_profile_component_time=0.0;
        double update_profile_uf2_time=0.0;
        double update_profile_connectivity_time=0.0;
        int update_profile_epsilon_num=0;
        int update_profile_ppt_num=0;
        long long update_profile_mu_num=0;
        long long update_profile_subgraph_vertex_num=0;
        long long update_profile_uf1_neighbor_num=0;
        long long update_profile_cs_vertex_num=0;
        long long update_profile_non_core_neighbor_num=0;
        long long update_profile_uf2_cs_vertex_num=0;
        long long update_profile_uf2_edge_num=0;
        long long update_profile_connectivity_neighbor_num=0;
        long long update_profile_new_edge_num=0;
        long long update_profile_replace_edge_num=0;
        long long update_profile_edge_num_before=0;
        long long update_profile_edge_num_after=0;
        vector<int> temp_label;
        UnionFind uf;
    public:
        SCAN_PPT_Index(string dataset, int DELTA=10, int format_id=-2)
        {
            map<int,string> formats = {
                {-1,"SortSet"},
                {-2,"Table"},
                {-3,"RStar"},
            };
            this->format=formats[format_id];
            this->DELTA=DELTA;
            cout<<"[INFORMATION] current ppt index's format: "<<format<<endl;
            init_graph(dataset);
            k=(int)(ceil((1/(2*pow(RHO,2)))*(log(2/FAILURE_PB))));
            this->uf=UnionFind(n);
        }
        ~SCAN_PPT_Index()
        {
            delete[] max_mu_under_epsilon;
        }
        void statistics_init()
        {
            traverse_node_num=0;
            traverse_core_num=0;
            traverse_non_core_num=0;
            traverse_edge_num=0;
            traverse_unite_num=0;
            traverse_pc_find_num=0;
            traverse_unite_time=0.0;
            traverse_pc_find_time=0.0;
        }
        void temp_label_init()
        {
            temp_label.assign(n, INF);
            uf.Init();
        }
        static bool cmp_rnk(const int& u, const int& v) 
        {
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
            interval_num.resize(n);
            for(int u=0;u<n;u++)
            {
                for(int j=0;j<DELTA;j++)
                {
                    interval_num[u].push_back(0);
                }
            }
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
                        int bucket=similarity/(PRECISION/DELTA);
                        if(bucket==DELTA) bucket=DELTA-1;
                        interval_num[u][bucket]+=1;
                        interval_num[v][bucket]+=1;
                    }
                }
                sort(neighbor_order[u].begin(),neighbor_order[u].end(),cmp_order_space_index);
            }
            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish building Neighbor Order"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            total_build_time+=construct_time.count();
        }
        void fast_construct_neighbor_order()
        {
            interval_num.resize(n);
            for(int u=0;u<n;u++)
            {
                for(int j=0;j<DELTA;j++)
                {
                    interval_num[u].push_back(0);
                }
            }
            string neighbor_order_file="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            ifstream read_index;
            read_index.open(neighbor_order_file);
            string line;
            int v_id=0;
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
                    int bucket=sim/(PRECISION/DELTA);
                    if(bucket==DELTA) bucket=DELTA-1;
                    interval_num[v_id][bucket]+=1;
                }
                neighbor_order.push_back(item);
                v_id++;
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
                    construct_time += t;
                    read_time=false;
                }
            }
            read_index.close();
            cout<<"[CONSTRUCT INDEX] finish loading Neighbor Order"<<endl;
            cout<<"[TIME COST] "<<construct_time<<" s."<<endl;
            total_build_time+=construct_time;
        }
        void construct_ppt_vertex_map()
        {
            cout<<"[CONSTRUCT Index] start PPT decomposition"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            vertex2ppt.resize(n);
            max_mu_under_epsilon=new int[DELTA];
            for(int i=0;i<DELTA;i++)
            {
                max_mu_under_epsilon[i]=1;
            }
            for(int i=0;i<n;i++)
            {
                int count=1;
                for(int bucket=DELTA-1;bucket>0;--bucket)
                {
                    count+=interval_num[i][bucket];
                    if(interval_num[i][bucket]!=0)
                    {
                        // cout<<" now vertex "<<i<<" belongs to ("<<bucket<<","<<count<<")-core"<<endl;
                        
                        pair<int,int> temp_point=make_pair(bucket,count);
                        auto it=point2ppt.find(temp_point);
                        //if PPT is not exist
                        if(it==point2ppt.end())
                        {
                            //construct PPT
                            PPT temp_PPT;
                            temp_PPT.id=ppt_max_id;  
                            temp_PPT.epsilon=bucket;
                            temp_PPT.mu=count;
                            temp_PPT.vertices.push_back(i);

                            //add PPT into PPT_set and assign a id, map PPT and id 
                            ppts.push_back(temp_PPT);
                            point2ppt.emplace(temp_point,ppt_max_id);
                                                        
                            //update max mu under each epsilon
                            if(max_mu_under_epsilon[bucket]<count)
                            {
                                max_mu_under_epsilon[bucket]=count;
                            }

                            //add PPT into vertex2ppt
                            vertex2ppt[i].push_back(ppt_max_id);     
                            
                            ppt_max_id++;
                        }
                        else
                        {
                            int PPT_id=it->second;

                            ppts[PPT_id].vertices.push_back(i);

                            vertex2ppt[i].push_back(PPT_id);
                        }
                    }
                }        
            }
            ppt_matrix.resize(DELTA);
            for(int i=0;i<DELTA;i++)
            {
                for(int j=2;j<max_mu_under_epsilon[i]+1;j++)
                {
                    ppt_matrix[i].push_back(-1);
                }
            }
            for(int i=0;i<ppts.size();i++)
            {
                ppt_matrix[ppts[i].epsilon][ppts[i].mu-2]=i;
            }
            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish PPT decomposition"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            total_build_time+=construct_time.count();
        }
        void construct_componets(int max_update_epsilon=-1)
        {
            cout<<"[CONSTRUCT INDEX] start building componets"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            long long profile_epsilon_num=0;
            long long profile_mu_num=0;
            long long profile_subgraph_vertex_num=0;
            long long profile_uf1_neighbor_num=0;
            long long profile_cs_vertex_num=0;
            long long profile_non_core_neighbor_num=0;
            long long profile_uf2_cs_vertex_num=0;
            long long profile_uf2_edge_num=0;
            long long profile_connectivity_neighbor_num=0;
            long long profile_new_edge_num=0;
            long long profile_replace_edge_num=0;
            if(max_update_epsilon<0)
            {
                max_update_epsilon=DELTA-1;
                cs2ppt.clear();
                vertex2cs.clear();
                vertex2cs.resize(n);
                for(int vertex=0;vertex<n;vertex++)
                {
                    vertex2cs[vertex].resize(vertex2ppt[vertex].size(),-1);
                }
            }
            //init
            vector<int> pre_subgraph;
            pre_subgraph.reserve(n);
            vector<int> pre_non_cores;
            pre_non_cores.reserve(n);
            vector<bool> vertex_bit_set;
            vertex_bit_set.resize(n);
            vector<bool> non_core_bit_set;
            non_core_bit_set.resize(n);
            vector<int> vertex2currentPPT(n,-1);
            vector<int> vertex2currentcs(n,-1);
            std::fill_n(non_core_bit_set.begin(), non_core_bit_set.size(), false);
            std::fill_n(vertex_bit_set.begin(), vertex_bit_set.size(), false);
            UnionFind UF(n);            //the entire connectivity in dominate space
            UnionFind UF2(n);           //the connectivity in pre-level dominate space

            for(int temp_eps=max_update_epsilon;temp_eps>=0;--temp_eps)
            {
                //if doesn't exist
                if(!ppt_matrix[temp_eps].size()) continue;
                profile_epsilon_num++;
                // cout<<"================================================================================="<<endl;
                // cout<<"now temp_eps: "<<temp_eps<<endl;

                //obtain connectivity range
                int temp_range=temp_eps*(PRECISION/DELTA);

                //init 
                UF.Init(pre_subgraph);
                UF2.Init(pre_subgraph);
                for(auto & v:pre_subgraph)
                {
                    vertex_bit_set[v]=false;
                    vertex2currentcs[v]=-1;
                    vertex2currentPPT[v]=-1;
                }
                for(auto & v:pre_non_cores)
                {non_core_bit_set[v]=false;}
                pre_subgraph.clear();
                pre_non_cores.clear();
                
                // fix temp_eps, PPT(temp_eps,temp_mu+1) must dominate PPT(temp_eps,temp_mu);
                for(int temp_mu=max_mu_under_epsilon[temp_eps];temp_mu>1;--temp_mu)
                {
                    profile_mu_num++;
                    //if not construct mcr-tree
                    //if not such PPT (temp_eps,temp_mu)
                    // if(ppt_matrix[temp_eps][temp_mu-2]==-1) continue;

                    // cout<<"now temp_mu: "<<temp_mu<<endl;
                    vector<int> subgraph;
                    vector<int> dominate_space;
                        
                    //find dominate space and its vertices, time complexity: O(DELTA)
                    int current_ppt_id=ppt_matrix[temp_eps][temp_mu-2];
                    if(current_ppt_id!=-1) dominate_space.push_back(current_ppt_id);

                    for(int temp_larger_eps=temp_eps+1;temp_larger_eps<DELTA;temp_larger_eps++)
                    {
                        if(max_mu_under_epsilon[temp_larger_eps]>=temp_mu)
                        {
                            int temp_ppt_id=ppt_matrix[temp_larger_eps][temp_mu-2];
                            if(temp_ppt_id!=-1)
                            {dominate_space.push_back(temp_ppt_id);}
                        }
                    }

                    //if current PPT is (temp_eps,max_mu_under_temp_eps);
                    if(temp_mu==max_mu_under_epsilon[temp_eps])
                    {
                        //add all ppt(>temp_eps,>max_mu_under_temp_eps) into dominate space
                        for(int temp_larger_eps=temp_eps+1;temp_larger_eps<DELTA;temp_larger_eps++)
                        {
                            if(max_mu_under_epsilon[temp_larger_eps]>temp_mu)
                            {
                                for(int temp_larger_mu=temp_mu+1;temp_larger_mu<=max_mu_under_epsilon[temp_larger_eps];temp_larger_mu++)
                                {
                                    int temp_ppt_id=ppt_matrix[temp_larger_eps][temp_larger_mu-2];
                                    if(temp_ppt_id!=-1)
                                    {dominate_space.push_back(temp_ppt_id);}
                                }
                            }
                        }
                    }

                    //add vertice to pre_subgraph and subgraph
                    for(auto sp_id:dominate_space)
                    {
                        for(auto v:ppts[sp_id].vertices)
                        {
                            if(!vertex_bit_set[v])
                            {
                                pre_subgraph.push_back(v);
                                subgraph.push_back(v);
                                vertex2currentPPT[v]=sp_id;
                                vertex_bit_set[v]=true;
                            }
                        }
                        for(auto cs_id : ppts[sp_id].componet_core_ids)
                        {
                            for(auto v: non_cores[cs_id])
                            {
                                if(!non_core_bit_set[v])
                                {
                                    pre_non_cores.push_back(v);
                                    non_core_bit_set[v]=true;
                                }
                            }
                        }
                    }
                    profile_subgraph_vertex_num+=subgraph.size();

                    // cout<<"dominate space: ";
                    // for(auto sp_id:dominate_space)
                    // {
                    //     cout<<"("<<ppts[sp_id].epsilon<<","<<ppts[sp_id].mu<<") ";
                    // }
                    // cout<<endl;
                    // cout<<"pre_subgraph: ";
                    // for(auto v:pre_subgraph)
                    // {cout<<v<<" ";}
                    // cout<<endl;
                    // cout<<"subgraph: ";
                    // for(auto v:subgraph)
                    // {cout<<v<<" ";}
                    // cout<<endl;
                    // cout<<"pre_non_cores: ";
                    // for(auto v:pre_non_cores)
                    // {
                    //     cout<<v<<' ';
                    // }
                    // cout<<endl;

                    //compute connectivity with union-find disjoint set
                    for(auto & v:subgraph)
                    {
                        for(auto & nv:neighbor_order[v])
                        {
                            //sigma(v,nv)>=range and nv is in now subgraph and root(v)!=root(nv);
                            if(nv.first<temp_range) break;
                            profile_uf1_neighbor_num++;
                            if(vertex_bit_set[nv.second]&&UF.Find(v)!=UF.Find(nv.second))
                            {UF.Unite(v,nv.second);}
                        }
                    }

                    //divide componet for this skyline parameter pair
                    if(current_ppt_id!=-1)
                    {
                        unordered_map<int,vector<int>> root2vertex;
                        profile_cs_vertex_num+=ppts[current_ppt_id].vertices.size();
                        for(auto & v:ppts[current_ppt_id].vertices)
                        {root2vertex[UF.Find(v)].push_back(v);}
                        for(auto & item:root2vertex)
                        {
                            int cs_id=cluster_slices.size();
                            ppts[current_ppt_id].componet_core_ids.push_back(cs_id);
                            cluster_slices.push_back(item.second);
                            cs2ppt.push_back(current_ppt_id);

                            for(auto vertex:item.second)
                            {
                                if(vertex2cs[vertex].size()<vertex2ppt[vertex].size())
                                {
                                    vertex2cs[vertex].resize(vertex2ppt[vertex].size(),-1);
                                }
                                for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                                {
                                    if(vertex2ppt[vertex][i]==current_ppt_id)
                                    {
                                        vertex2cs[vertex][i]=cs_id;
                                        break;
                                    }
                                }
                            }
                        }
                        for(auto cs_id:ppts[current_ppt_id].componet_core_ids)
                        {
                            vector<int> temp_non_cores;
                            // unordered_set<int> record;
                            for(auto v:cluster_slices[cs_id])
                            {
                                for(auto & nv:neighbor_order[v])
                                {
                                    if(nv.first<temp_range) break;
                                    profile_non_core_neighbor_num++;
                                    if(!vertex_bit_set[nv.second]&&!non_core_bit_set[nv.second])
                                    {
                                        // record.insert(nv.second);
                                        temp_non_cores.push_back(nv.second);
                                        pre_non_cores.push_back(nv.second);
                                        non_core_bit_set[nv.second]=true;
                                    }
                                }
                            }
                            // vector<int> temp_non_cores(record.begin(),record.end());
                            non_cores.push_back(temp_non_cores);
                        }
                    }
                    
                    //connect the componets that has connection
                    for(auto & sp_id:dominate_space)
                    {
                        for(auto cs_id:ppts[sp_id].componet_core_ids)
                        {
                            profile_uf2_cs_vertex_num+=cluster_slices[cs_id].size();
                            auto it = cluster_slices[cs_id].begin();
                            int firstVertex = *it;
                            vertex2currentcs[*it]=cs_id;
                            it++;
                            while (it != cluster_slices[cs_id].end()) {
                                UF2.Unite(firstVertex, *it);
                                vertex2currentcs[*it]=cs_id;
                                it++;
                            }
                        }
                        profile_uf2_edge_num+=ppts[sp_id].edges.size();
                        for(auto & edge:ppts[sp_id].edges)
                        {UF2.Unite(edge.first, edge.second);}
                    }
                    //compute the new connection between componets
                    // if(current_ppt_id!=-1)
                    // {
                    for(auto & v:subgraph)
                    {
                        // cout<<"now deal "<<v<<endl;
                        auto v_ppt_id=vertex2currentPPT[v];
                        // cout<<v<<" corresponding to PPT ("<<ppts[v_ppt_id].epsilon<<","<<ppts[v_ppt_id].mu<<")"<<endl;
                        for(auto & nv:neighbor_order[v])
                        {
                            if(nv.first<temp_range) break;
                            profile_connectivity_neighbor_num++;
                            // if(vertex2currentPPT.find(nv.second)==vertex2currentPPT.end()) continue;
                            if(vertex_bit_set[nv.second]&&UF2.Find(v)!=UF2.Find(nv.second))
                            {
                                // cout<<"now deal its neighbor "<<nv.second<<endl;
                                auto temp_cs_id=vertex2currentcs[nv.second];
                                UF2.Unite(v,nv.second);
                                // if e is the connection between current PPT and nv's PPT, we store the edge of max{sigma(e),sigma(v,nv)}
                                auto it=ppts[v_ppt_id].cs2edges.find(temp_cs_id);
                                if(it!=ppts[v_ppt_id].cs2edges.end())
                                {
                                    int edge_id=it->second;
                                    if(nv.first>ppts[v_ppt_id].connectivitys[edge_id])
                                    {
                                        profile_replace_edge_num++;
                                        ppts[v_ppt_id].edges[edge_id]=make_pair(v,nv.second);
                                        ppts[v_ppt_id].connectivitys[edge_id]=nv.first;
                                    }
                                }
                                else
                                {
                                    profile_new_edge_num++;
                                    ppts[v_ppt_id].cs2edges[temp_cs_id]=ppts[v_ppt_id].edges.size();
                                    ppts[v_ppt_id].edges.push_back(make_pair(v,nv.second));
                                    ppts[v_ppt_id].connectivitys.push_back(nv.first);
                                }
                            }
                        }
                    }
                }
            }
            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish building componets"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            cout<<"[CONSTRUCT PROFILE] epsilon: "<<profile_epsilon_num
                <<" mu: "<<profile_mu_num
                <<" subgraph vertices: "<<profile_subgraph_vertex_num
                <<" UF1 neighbor visits: "<<profile_uf1_neighbor_num<<endl;
            cout<<"[CONSTRUCT PROFILE] CS vertices: "<<profile_cs_vertex_num
                <<" non-core neighbor visits: "<<profile_non_core_neighbor_num
                <<" UF2 CS vertices: "<<profile_uf2_cs_vertex_num
                <<" UF2 edges: "<<profile_uf2_edge_num<<endl;
            cout<<"[CONSTRUCT PROFILE] connectivity neighbor visits: "
                <<profile_connectivity_neighbor_num
                <<" new edges: "<<profile_new_edge_num
                <<" replaced edges: "<<profile_replace_edge_num<<endl;
            total_build_time+=construct_time.count();
            for(auto & ppt:ppts)
            {
                total_build_space+=12;
                total_build_space+=ppt.componet_core_ids.size()*4;
                edge_space+=ppt.edges.size()*8;
                connectivity_space+=ppt.connectivitys.size()*4;
                total_build_space+=ppt.edges.size()*12;
            }
            for(auto & item:cluster_slices)
            {
                total_build_space+=item.size()*4;
            }
            for(auto & item:non_cores)
            {
                total_build_space+=item.size()*4;
            }
        }
        void construct_kd_tree()
        {

        }
        void construct_r_star_tree()
        {

        }
        void construct() override
        {
            #ifdef FAST_CONSTRUCT
                fast_construct_neighbor_order();
            #else
                construct_bottom_k_sketch();
                construct_neighbor_order();
            #endif
            construct_ppt_vertex_map();
            construct_componets();
        }
        void connect_vertices_after_search_index(int bucket, int mu, vector<int> &ppt_ids) 
        {
            // const auto begin=std::chrono::steady_clock().now();  
            vector<int> all_cs_ids;
            for (int ppt_id : ppt_ids) 
            {
                for (int cs_id : ppts[ppt_id].componet_core_ids)
                {
                    all_cs_ids.push_back(cs_id);
                }
            }
            for (auto cs_id : all_cs_ids) {
                if (cluster_slices[cs_id].empty()) continue;
                int rep = *cluster_slices[cs_id].begin();   
                for (int v : cluster_slices[cs_id]) {
                    temp_label[v] = rep;  
                }
            }
            const auto end1=std::chrono::steady_clock().now();
            for (auto ppt_id : ppt_ids) {
                for (int i = 0; i < (int)ppts[ppt_id].edges.size(); i++) 
                {
                    if (ppts[ppt_id].connectivitys[i] >= bucket*(PRECISION/DELTA)) 
                    {
                        auto &e = ppts[ppt_id].edges[i];
                        int rep_u = temp_label[e.first];
                        int rep_v = temp_label[e.second];
                        if(rep_u == INF || rep_v == INF || rep_u == rep_v) continue;
                        // traverse_edge_num++;
                        // auto unite_begin=std::chrono::steady_clock().now();
                        uf.Unite(rep_u, rep_v);
                        // auto unite_end=std::chrono::steady_clock().now();
                        // tms unite_time=unite_end-unite_begin;
                        // traverse_unite_time+=unite_time.count();
                        // traverse_unite_num++;
                    }
                }
            }
            const auto end2=std::chrono::steady_clock().now();

            // vector<int> cluster_ids;
            // for (auto cs_id : all_cs_ids) 
            // {
            //     int rep_v = *cluster_slices[cs_id].begin();
            //     cluster_ids.push_back(uf.Find(rep_v));
            // }

            // int cluster_idx = 0;
            for (auto cs_id : all_cs_ids) {
                // int cluster_id = cluster_ids[cluster_idx++];
                int cluster_id=uf.Find(cluster_slices[cs_id][0]);

                for (auto core_v : cluster_slices[cs_id]) {
                    // traverse_core_num++;
                    if (labels[core_v]!=cluster_id)
                        labels[core_v]=cluster_id;
                }
                for (auto non_core_v : non_cores[cs_id]) {
                    // traverse_non_core_num++;
                    if (labels[non_core_v] == INF)
                        labels[non_core_v] = cluster_id;
                }
            }
            // const auto end3=std::chrono::steady_clock().now();
            // tms stage1_time=end1-begin;
            // tms stage2_time=end2-end1;
            // tms stage3_time=end3-end2;
            // cout<<"stage1 time: "<<stage1_time.count()<<" stage2 time: "<<stage2_time.count()<<" stage3 time: "<<stage3_time.count()<<endl;
        }
        void query_from_sortset(int bucket,int mu)
        {
            cout<<"[CLUSTER] start query clusters from sorted set"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            vector<int> ppt_ids;
            ppt_ids.reserve(ppts.size());
            for(auto &ppt:ppts)
            {
                if(ppt.epsilon<bucket) break;
                if(ppt.mu>=mu)
                {
                    ppt_ids.push_back(ppt.id);
                }
            }
            #ifdef TEST_CORRECTNESS
                cout<<"[TEST] print collecting ppt: "<<endl;
                for(auto id:ppt_ids)
                {
                    cout<<"ppt id: "<<id<<" epsilon: "<<ppts[id].epsilon<<" mu: "<<ppts[id].mu<<endl;
                }
            #endif
            #ifdef TEST_STRUCTURE
                const auto end=std::chrono::steady_clock().now();
                tms collection_time=end-begin;
                cout<<"[TEST] finish collecting PPT"<<endl;
                cout<<"[TIME COST] "<<collection_time.count()<<" s."<<endl;
                total_collection_ppt_time+=collection_time.count();
            #else
                connect_vertices_after_search_index(bucket,mu,ppt_ids);
                const auto end=std::chrono::steady_clock().now();
                tms cluster_time=end-begin;
                cout<<"[CLUSTER] finish clustering"<<endl;
                cout<<"[TIME COST] "<<cluster_time.count()<<" s."<<endl;
                total_cluster_time+=cluster_time.count();
            #endif
        }
        void query_from_table(int bucket,int mu)
        {
            cout<<"[CLUSTER] start query clusters from table"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            vector<int> ppt_ids;
            ppt_ids.reserve(ppts.size());
            for(int i=table_ppts.size()-1;i>=bucket;--i)
            {
                const auto& vec = table_ppts[i];
                for (int ppt_id : vec)
                {
                    if (ppts[ppt_id].mu < mu) break;
                    ppt_ids.push_back(ppt_id);
                }
            }
            #ifdef TEST_CORRECTNESS
                cout<<"[TEST] print collecting ppt: "<<endl;
                for(auto id:ppt_ids)
                {
                    cout<<"ppt id: "<<id<<" epsilon: "<<ppts[id].epsilon<<" mu: "<<ppts[id].mu<<endl;
                }
            #endif
            #ifdef TEST_STRUCTURE
                const auto end=std::chrono::steady_clock().now();
                tms collection_time=end-begin;
                cout<<"[TEST] finish collecting PPT"<<endl;
                cout<<"[TIME COST] "<<collection_time.count()<<" s."<<endl;
                total_collection_ppt_time+=collection_time.count();
            #else
                connect_vertices_after_search_index(bucket,mu,ppt_ids);
                const auto end=std::chrono::steady_clock().now();
                tms cluster_time=end-begin;
                cout<<"[CLUSTER] finish clustering"<<endl;
                cout<<"[TIME COST] "<<cluster_time.count()<<" s."<<endl;
                total_cluster_time+=cluster_time.count();
            #endif
        }
        void query_from_r_star_tree(int bucket,int mu)
        {
            cout<<"[CLUSTER] start query clusters from table"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            vector<int> ppt_ids;
            ppt_ids.reserve(ppts.size());
            
            auto bound=bounds(bucket,mu,DELTA,dmax);
            r_star_tree.Query(RTree::AcceptEnclosing(bound),Visitor(),ppt_ids);
            #ifdef TEST_CORRECTNESS
                cout<<"[TEST] print collecting ppt: "<<endl;
                for(auto id:ppt_ids)
                {
                    cout<<"ppt id: "<<id<<" epsilon: "<<ppts[id].epsilon<<" mu: "<<ppts[id].mu<<endl;
                }
            #endif
            #ifdef TEST_STRUCTURE
                const auto end=std::chrono::steady_clock().now();
                tms collection_time=end-begin;
                cout<<"[TEST] finish collecting PPT"<<endl;
                cout<<"[TIME COST] "<<collection_time.count()<<" s."<<endl;
                total_collection_ppt_time+=collection_time.count();
            #else
                connect_vertices_after_search_index(bucket,mu,ppt_ids);
                const auto end=std::chrono::steady_clock().now();
                tms cluster_time=end-begin;
                cout<<"[CLUSTER] finish clustering"<<endl;
                cout<<"[TIME COST] "<<cluster_time.count()<<" s."<<endl;
                total_cluster_time+=cluster_time.count();
            #endif
        }
        int locate(int epsilon)
        {
            int bucket_range=PRECISION/DELTA;
            int bucket_id= epsilon/bucket_range;
            return bucket_id;
        }
        void check_query(int bucket, int mu)
        {
            cout<<"[DEBUG] query from PPT Index "<<endl;
            vector<int> ppt_ids;
            ppt_ids.reserve(ppts.size());
            for(int i=table_ppts.size()-1;i>0;--i)
            {
                if(bucket>i) break;
                for(int j=0;j<table_ppts[i].size();j++)
                {
                    traverse_node_num++;
                    int ppt_id=table_ppts[i][j];
                    if(ppts[ppt_id].mu<mu) break;
                    ppt_ids.push_back(ppt_id);
                }
            }
            cout<<"containing "<<ppt_ids.size()<<" PPTS: "<<endl;
            vector<int> core_labels(n,-1);
            unordered_map<int,vector<int>> clusters;
                         
            for (auto ppt_id : ppt_ids) {
                for (auto cs_id : ppts[ppt_id].componet_core_ids) {
                    if (cluster_slices[cs_id].empty()) continue;
                    int rep = *cluster_slices[cs_id].begin();   
                    for (int v : cluster_slices[cs_id]) {
                        if(temp_label[v]!=-1)
                        {uf.Unite(temp_label[v],rep);}
                        else{temp_label[v] = rep;}
                    }
                }
            }
            for (auto ppt_id : ppt_ids) {
                for (int i = 0; i < (int)ppts[ppt_id].edges.size(); i++) 
                {
                    if (ppts[ppt_id].connectivitys[i] >= bucket*(PRECISION/DELTA)) 
                    {
                        auto &e = ppts[ppt_id].edges[i];
                        int rep_u = temp_label[e.first];
                        int rep_v = temp_label[e.second];
                        if(rep_u == INF || rep_v == INF || rep_u == rep_v) continue;
                        uf.Unite(rep_u, rep_v);
                    }
                }
            }
            vector<int> cluster_ids;
            for (auto ppt_id : ppt_ids) {
                for (auto cs_id:ppts[ppt_id].componet_core_ids) {
                    int rep_v = *cluster_slices[cs_id].begin();
                    cluster_ids.push_back(uf.Find(rep_v));
                }
            }

            int cluster_idx = 0;
            for (auto &ppt_id : ppt_ids) {
                for (auto cs_id:ppts[ppt_id].componet_core_ids) {
                    int cluster_id = cluster_ids[cluster_idx++];
                    if(clusters.find(cluster_id)==clusters.end())
                    {
                        clusters[cluster_id]=cluster_slices[cs_id];
                    }
                    else
                    {
                        for (auto core_v : cluster_slices[cs_id]) 
                        {
                            clusters[cluster_id].push_back(core_v);
                        }
                    }
                }
            }
            for(auto & item:clusters)
            {
                int min_id=n;
                for(auto v:item.second)
                {
                    if(v<min_id) min_id=v;
                }
                for(auto v:item.second)
                {
                    core_labels[v]=min_id;
                }
            }
            cout<<"[DEBUG] clustering of PPT finish, cluster num: "<<clusters.size()<<endl;
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

            vector<vector<pair<int, int>>> bucket_index;
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
                bucket_index.push_back(item);
            }
            read_index.close();
            cout<<"[LOAD Index] finish loading no and ci !!!!"<<endl;
            cout<<"[CLUSTER] clustering with BOTBIN start !!!"<<endl;
            cout<<"search from bucket["<<bucket<<"]"<<endl;
            vector<int> core_label_truth(n,-1);
            unordered_set<int> core_vertice_set;
            for(auto &core:bucket_index[bucket])
            {
                if(core.first<mu-1) break;
                core_vertice_set.insert(core.second);
            }
            int cluster_num_of_botbin=0;
            for(auto core:bucket_index[bucket])
            {
                if(core_label_truth[core.second]!=-1) continue;
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
                        // cout<<"("<<nv.first<<","<<nv.second<<")"<<' ';
                        if(nv.first<bucket*(PRECISION/DELTA)) break;
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
                    core_label_truth[item]=min_id_cluster;
                }
                // cout<<endl;
            }
            cout<<"[CLSUTER] finish clustering of botbin, cluster results num: "<<cluster_num_of_botbin<<endl;
            for(int i=0;i<n;i++)
            {
                if(core_label_truth[i]==-1&&core_labels[i]!=-1)
                {
                    cout<<"vertex "<<i<<" is assign in a cluster of ppt but not in a cluster of botbin"<<endl;
                }
                else if(core_label_truth[i]!=-1&&core_labels[i]==-1)
                {
                    cout<<"vertex "<<i<<" is assign in a cluster of botbin but not in a cluster of ppt"<<endl;
                }
                else if(core_label_truth[i]==-1&&core_labels[i]==-1)
                {
                    continue;
                }
                else if(core_label_truth[i]!=core_labels[i])
                {
                    cout<<"vertex "<<i<<" belongs to v "<<core_label_truth[i]<<" of botbin and belongs to v "<<core_labels[i]<<" of ppt"<<endl;
                }
                else
                {
                    continue;
                }
            }
            cout<<"=================================================================="<<endl;
            cout<<"Print vertices that are used for debug "<<endl;
            unordered_set<int> watch_vertex={13987, 152106};
            for(auto v:watch_vertex)
            {
                cout<<"vertex "<<v<<" neighbors: ";
                for(auto nv:neighbor_order[v])
                {cout<<"("<<nv.first<<","<<nv.second<<")"<<' ';}
                cout<<endl;
            }
            for(auto ppt_id:ppt_ids)
            {
                for(auto cs_id:ppts[ppt_id].componet_core_ids)
                {
                    bool print=false;
                    for(auto v:cluster_slices[cs_id])
                    {
                        if(watch_vertex.find(v)!=watch_vertex.end())
                        {
                            print=true;
                            break;
                        }
                    }
                    if(print)
                    {
                        cout<<"epsilon: "<<ppts[ppt_id].epsilon<<" mu: "<<ppts[ppt_id].mu<<endl;
                        cout<<"edge and connectivty: "<<endl;
                        cout<<"connectivity: ";
                        for(int i=0;i<ppts[ppt_id].edges.size();i++)
                        {
                            if(watch_vertex.find(ppts[ppt_id].edges[i].first)!=watch_vertex.end()||watch_vertex.find(ppts[ppt_id].edges[i].second)!=watch_vertex.end())
                            {cout<<"("<<ppts[ppt_id].edges[i].first<<","<<ppts[ppt_id].edges[i].second<<"):"<<ppts[ppt_id].connectivitys[i]<<" ";}
                        }
                        cout<<endl;
                        cout<<"cluster slices "<<cs_id<<":";
                        for(auto v:cluster_slices[cs_id])
                        {
                            cout<<v<<' ';
                        }
                        cout<<endl;
                    }
                }
            }
        }
        void query(float epsilon, int mu) override
        {
            #ifdef CHECK_RESULT
                clear_cluster_result();
                statistics_init();
                temp_label_init();
                cout<<"[QUERY] now query epsilon: "<<epsilon<<" mu: "<<mu<<endl;
                int eps=epsilon*PRECISION;
                int bucket=locate(eps);
                check_query(bucket,mu);
            #else
                clear_cluster_result();
                statistics_init();
                temp_label_init();
                cout<<"[QUERY] now query epsilon: "<<epsilon<<" mu: "<<mu<<endl;
                int eps=epsilon*PRECISION;
                int bucket=locate(eps);
                if(format=="SortSet")
                {
                    query_from_sortset(bucket,mu);
                }
                if(format=="Table")
                {
                    query_from_table(bucket,mu);
                }
                else if(format=="KD")
                {

                }
                else if(format=="RStar")
                {
                    query_from_r_star_tree(bucket,mu);
                }
                else
                {

                }
                query_count++;
                #ifdef EVALUATION_CLUSTERING_QUALITY
                    evaluate_clustering_quality(epsilon,mu);
                #endif
            #endif
            // print_statistics();                                                                
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
            if(!read_groundtruth) cout<<"error! no such ground truth!!!"<<endl;
            vector<int> groundtruth(n);
            read_groundtruth.read(reinterpret_cast<char*>(groundtruth.data()),sizeof(int)*n);
            read_groundtruth.close();
            double temp_ari=adjusted_rand_index_ignore_minus1(groundtruth,labels);
            cout<<"[EVALUATION] ari index score: "<<temp_ari<<endl;
            total_ari+=temp_ari;
        }
        void print_all_PPT()
        {
            cout<<"[RPINT] now print all ppt after construct componets"<<endl;
            for(auto & ppt:ppts)
            {
                cout<<"================================================================================="<<endl;
                cout<<"id: "<<ppt.id<<" epsilon: "<<ppt.epsilon<<" mu: "<<ppt.mu<<endl;
                cout<<"vertices: ";
                for(auto & v:ppt.vertices)
                {
                    cout<<v<<' ';
                }
                cout<<endl;
                cout<<"subgraph: ";
                for(auto cs_id:ppt.componet_core_ids)
                {
                    cout<<"(";
                    cout<<"cores: ";
                    for(auto & v:cluster_slices[cs_id])
                    {cout<<v<<" ";}
                    cout<<"non-cores: ";
                    for(auto & v:non_cores[cs_id])
                    {cout<<v<<" ";}
                    cout<<")";
                }
                cout<<endl;
                cout<<"connectivity: ";
                for(int i=0;i<ppt.edges.size();i++)
                {
                    cout<<"("<<ppt.edges[i].first<<","<<ppt.edges[i].second<<"):"<<ppt.connectivitys[i]<<" ";
                }
                cout<<endl;
            }
        }
        void print_similarity_interval()
        {
            cout<<"[PRINT] now print the similarity interval"<<endl;
            for(int i=0;i<DELTA;i++)
            {
                cout<<"["<<float(i/float(DELTA))<<","<<float((i+1)/float(DELTA))<<") ";
            }
            cout<<endl;
            cout<<"[PRINT] print finish !!!"<<endl;
        }
        void print_PPT_vertex_map()
        {
            cout<<"[PRINT] now print the map between PPT and vertex !!!"<<endl;
            
            cout<<"All PPT: "<<endl;
            for(auto & PPT:ppts)
            {
                cout<<"("<<PPT.epsilon<<","<<PPT.mu<<")"<<" vertices: ";
                for(auto & v:PPT.vertices)
                {
                    cout<<v<<" ";
                }
                cout<<endl;
            }

            cout<<"max mu under each epsilon: "<<endl;
            for(int i=0;i<DELTA;i++)
            {
                cout<<"epsilon bucket "<<i<<" : "<<max_mu_under_epsilon[i]<<endl;
            }

            cout<<"PPT matrix: "<<endl;
            for(auto & line:ppt_matrix)
            {
                for(auto & elem:line)
                {
                    cout<<elem<<' ';
                }
                cout<<endl;
            }

            cout<<"Vertex 2 PPT: "<<endl;
            for(int i=0;i<vertex2ppt.size();i++)
            {
                cout<<"vertex: "<<i<<" PPT: ";
                for(auto & PPT_id:vertex2ppt[i])
                {
                    cout<<"("<<ppts[PPT_id].epsilon<<","<<ppts[PPT_id].mu<<")";
                }
                cout<<endl;
            }

            cout<<"[PRINT] print finish !!!"<<endl;
        }
        void print_statistics()
        {
            cout<<"[INFORMATION DEBUG] traverse node num: "<<traverse_node_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse core num: "<<traverse_core_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse non core num: "<<traverse_non_core_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse edge num: "<<traverse_edge_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse unite num: "<<traverse_unite_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse unite time: "<<traverse_unite_time<<endl;
            cout<<"[INFORMATION DEBUG] traverse path-compress find num: "<<traverse_pc_find_num<<endl;
            cout<<"[INFORMATION DEBUG] traverse path-compress find time: "<<traverse_pc_find_time<<endl;
        }
        void print_build_cost() override
        {
            cout<<"[PRINT NODE NUM] total point (leaf node): "<<ppts.size()<<endl;
            int total_core_componets=cluster_slices.size();
            double avg_core_componets=(double)total_core_componets/(double)ppts.size();
            cout<<"[PRINT CLUSTER SLICE]: total core componets: "<<total_core_componets<<endl;
            cout<<"[PRINT CLUSTER SLICE]: average core componets per PPT: "<<avg_core_componets<<endl;
            cout<<"[PRINT TIME] total building time: "<<total_build_time<<endl;
            cout<<"[PRINT SPACE INFORMATION] core vertices space cost: "<<(core_vertices_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] non core vertices space cost: "<<(non_core_vertices_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] edges space cost: "<<(edge_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] connectivity space cost: "<<(connectivity_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT PPT SPACE]: total space cost: "<<(total_build_space/1024)/1024<<".MB"<<endl;
        }
        void print_cluster_time() override
        {
            cout<<"[PRINT TIME] total clustering time: "<<total_cluster_time<<endl;
            cout<<"[PRINT TIME] average clustering time: "<<total_cluster_time/double(query_count)<<endl;
            #ifdef EVALUATION_CLUSTERING_QUALITY
                cout<<"[PRINT ARI] average ARI score: "<<total_ari/double(query_count)<<endl;
            #endif
            #ifdef TEST_STRUCTURE
                cout<<"[PRINT TIME] total collection time: "<<total_collection_ppt_time<<endl;
                cout<<"[PRINT TIME] average collection time: "<<total_collection_ppt_time/double(query_count)<<endl;
            #endif
        }
        void print_neighbor_order()
        {
            cout<<"[PRINT] print NEIGHBOR ORDER start !!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<i<<": ";
                for(auto & item:neighbor_order[i])
                {
                    cout<<"("<<item.first<<","<<item.second<<")";
                }
                cout<<endl;
            }
            cout<<"[PRINT] print NEIGHBOR ORDER finish !!!"<<endl;
        }
        void print_interval_num()
        {
            cout<<"[PRINT] print INTERVAL INDEX start !!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<i<<": ";
                for(auto & item:interval_num[i])
                {
                    cout<<item<<' ';
                }
                cout<<endl;
            }
            cout<<"[PRINT] print INTERVAL INDEX finish !!!"<<endl;
        }
        void print_index() override
        {
            // print_neighbor_order();
            print_similarity_interval();
            // print_interval_num();
            // print_PPT_vertex_map();
            print_all_PPT();
        }
        void load_index() override
        {
            cout<<"[LOAD INDEX] start !!!!"<<endl;
            string store_ppts="./Index/"+dataset+"/ppts_"+to_string(DELTA)+".txt";
            cout<<"read path: "<<store_ppts<<endl;
            ifstream read_index;
            string line;

            read_index.open(store_ppts);
            int ppt_size;
            int cs_size;
            read_index>>ppt_size>>cs_size;
            ppts.clear();
            ppts.resize(ppt_size);
            cluster_slices.clear();
            cluster_slices.resize(cs_size);
            non_cores.clear();
            non_cores.resize(cs_size);
            //read ppts
            for(int i=0;i<ppt_size;i++)
            {
                PPT & ppt=ppts[i];
                int componet_size;
                int edge_size;
                read_index>>ppt.id>>ppt.epsilon>>ppt.mu>>componet_size>>edge_size;
                ppt.componet_core_ids.resize(componet_size);
                ppt.edges.resize(edge_size);
                ppt.connectivitys.resize(edge_size);
                for(int j=0;j<componet_size;j++)
                {
                    read_index>>ppt.componet_core_ids[j];
                }
                for(int j=0;j<edge_size;j++)
                {
                    read_index>>ppt.edges[j].first>>ppt.edges[j].second;
                }
                for(int j=0;j<edge_size;j++)
                {
                    read_index>>ppt.connectivitys[j];
                }
            }
            //read cs and noncores
            for(int i=0;i<cs_size;i++)
            {
                int temp_size;
                auto & temp_item=cluster_slices[i];
                read_index>>temp_size;
                temp_item.resize(temp_size);
                for(int j=0;j<temp_size;j++)
                {
                    read_index>>temp_item[j];
                }
            }
            for(int i=0;i<cs_size;i++)
            {
                int temp_size;
                auto & temp_item=non_cores[i];
                read_index>>temp_size;
                temp_item.resize(temp_size);
                for(int j=0;j<temp_size;j++)
                {
                    read_index>>temp_item[j];
                }
            }
            read_index.close();
            if(format=="SortSet")
            {
                sort(ppts.begin(),ppts.end(),cmp_ppt);
                for(int i=0;i<ppts.size();i++)
                {
                    ppts[i].id=i;
                }
            }
            else if(format=="Table")
            {
                table_ppts.resize(DELTA);
                sort(ppts.begin(),ppts.end(),cmp_ppt);
                for(int i=0;i<ppts.size();i++)
                {
                    ppts[i].id=i;
                }
                for(auto ppt:ppts)
                {
                    table_ppts[ppt.epsilon].push_back(ppt.id);
                }
                for(int i=0;i<DELTA;i++)
                {
                    if(table_ppts[i].size())
                    {
                        sort(table_ppts[i].begin(),table_ppts[i].end(),[this](int a,int b){
                            return ppts[a].mu > ppts[b].mu;
                    });
                    }
                }
                int table_num=0;
                for(auto t:table_ppts)
                {
                    table_num+=t.size();
                }
                
            }
            else if(format=="KD")
            {

            }
            else if(format=="RStar")
            {
                sort(ppts.begin(),ppts.end(),cmp_ppt);
                for(int i=0;i<ppts.size();i++)
                {
                    ppts[i].id=i;
                }
                for(auto &ppt:ppts)
                {
                    r_star_tree.Insert(ppt.id,bounds(ppt.epsilon,ppt.mu,0,0));
                }
            }
            else{}
            // cout<<"================================================="<<endl;
            // print_all_PPT();
            int ppt_num=ppts.size();
            restore_ppt_membership_for_update();
            cout<<"ppt num: "<<ppt_num<<endl;
            cout<<"[LOAD INDEX] finish !!!!"<<endl;
        }
        void store_index() override
        {
            cout<<"[STORE INDEX] start !!!!"<<endl;
            string store_ppts="./Index/"+dataset+"/ppts_"+to_string(DELTA)+".txt";
            ofstream write_index;
            write_index.open(store_ppts);
            write_index<<ppts.size()<<" "<<cluster_slices.size()<<endl;
            for(auto & ppt:ppts)
            {
                write_index<<ppt.id<<' '<<ppt.epsilon<<' '<<ppt.mu<<' ';
                write_index<<ppt.componet_core_ids.size()<<' ';
                write_index<<ppt.edges.size()<<' ';
                for(auto cs_id:ppt.componet_core_ids)
                {
                    write_index<<cs_id<<' ';
                }
                for(auto edge:ppt.edges)
                {
                    write_index<<edge.first<<' '<<edge.second<<' ';
                }
                for(auto connectiviy:ppt.connectivitys)
                {
                    write_index<<connectiviy<<' ';
                }
                write_index<<endl;
            }
            for(auto & item:cluster_slices)
            {
                sort(item.begin(),item.end());
                write_index<<item.size()<<' ';
                for(auto v: item)
                {
                    write_index<<v<<' ';
                }
                write_index<<endl;
            }
            for(auto & item:non_cores)
            {
                sort(item.begin(),item.end());
                write_index<<item.size()<<' ';
                for(auto v: item)
                {
                    write_index<<v<<' ';
                }
                write_index<<endl;
            }
            write_index.close();
            cout<<"[STORE INDEX] finish !!!!"<<endl;
        }

        int compute_similarity_exact_for_update(int u, int v)
        {
            int common_neighbor_num=2;
            int i=0;
            int j=0;
            while(i<degree[u]&&j<degree[v])
            {
                if(graph[u][i]==graph[v][j])
                {
                    common_neighbor_num++;
                    i++;
                    j++;
                }
                else if(graph[u][i]>graph[v][j])
                {
                    i++;
                }
                else
                {
                    j++;
                }
            }

            double similarity=double(common_neighbor_num)/
                double(degree[u]+degree[v]+2-common_neighbor_num);
            return int(similarity*PRECISION);
        }

        void clear_ppt_structure_for_update()
        {
            ppts.clear();
            cluster_slices.clear();
            non_cores.clear();
            point2ppt.clear();
            vertex2ppt.clear();
            ppt_matrix.clear();
            table_ppts.clear();
            free_update_cs_ids.clear();
            ppt_max_id=0;

            delete[] max_mu_under_epsilon;
            max_mu_under_epsilon=nullptr;

            core_vertices_space=0.0;
            non_core_vertices_space=0.0;
            edge_space=0.0;
            connectivity_space=0.0;
            total_build_space=0.0;
        }

        void prepare_query_structure_after_update()
        {
            if(format=="SortSet")
            {
                sort(ppts.begin(),ppts.end(),cmp_ppt);
                for(int i=0;i<(int)ppts.size();i++)
                {
                    ppts[i].id=i;
                }
            }
            else if(format=="Table")
            {
                table_ppts.clear();
                table_ppts.resize(DELTA);
                for(int i=0;i<(int)ppts.size();i++)
                {
                    ppts[i].id=i;
                    table_ppts[ppts[i].epsilon].push_back(i);
                }
                for(int i=0;i<DELTA;i++)
                {
                    sort(table_ppts[i].begin(),table_ppts[i].end(),[this](int a,int b)
                    {
                        return ppts[a].mu>ppts[b].mu;
                    });
                }
            }
        }

        void update_similarity_graph_for_ppt(string update_type, int u, int v,
                                             vector<int> &affected_vertices,
                                             vector<vector<pair<int,int>>> &affected_mu_intervals)
        {
            int bucket_width=PRECISION/DELTA;
            unordered_set<int> dirty_neighbor_order;
            unordered_set<int> affected_vertex_set;
            vector<ppt_update_edge_t> changed_edges;
            dirty_neighbor_order.insert(u);
            dirty_neighbor_order.insert(v);
            affected_mu_intervals.clear();
            affected_mu_intervals.resize(DELTA);

            // The graph is already in its new state, while neighbor_order and
            // interval_num still describe the state before this edge update.
            int endpoints[2]={u,v};
            for(int k=0;k<2;k++)
            {
                int endpoint=endpoints[k];
                for(auto &item:neighbor_order[endpoint])
                {
                    int neighbor=item.second;
                    if(neighbor==u||neighbor==v) continue;

                    int old_similarity=item.first;
                    int new_similarity=compute_similarity_exact_for_update(endpoint,neighbor);
                    if(old_similarity==new_similarity) continue;

                    int old_bucket=old_similarity/bucket_width;
                    int new_bucket=new_similarity/bucket_width;
                    if(old_bucket==DELTA) old_bucket=DELTA-1;
                    if(new_bucket==DELTA) new_bucket=DELTA-1;
                    if(old_bucket!=new_bucket)
                    {
                        interval_num[endpoint][old_bucket]--;
                        interval_num[neighbor][old_bucket]--;
                        interval_num[endpoint][new_bucket]++;
                        interval_num[neighbor][new_bucket]++;
                        affected_vertex_set.insert(endpoint);
                        affected_vertex_set.insert(neighbor);
                        changed_edges.push_back({endpoint,neighbor,old_similarity,new_similarity});
                    }

                    item.first=new_similarity;
                    for(auto &reverse_item:neighbor_order[neighbor])
                    {
                        if(reverse_item.second==endpoint)
                        {
                            reverse_item.first=new_similarity;
                            break;
                        }
                    }
                    dirty_neighbor_order.insert(neighbor);
                }
            }

            if(update_type=="insert")
            {
                int similarity=compute_similarity_exact_for_update(u,v);
                int bucket=similarity/bucket_width;
                if(bucket==DELTA) bucket=DELTA-1;

                neighbor_order[u].push_back(make_pair(similarity,v));
                neighbor_order[v].push_back(make_pair(similarity,u));
                interval_num[u][bucket]++;
                interval_num[v][bucket]++;
                if(bucket>0)
                {
                    affected_vertex_set.insert(u);
                    affected_vertex_set.insert(v);
                    changed_edges.push_back({u,v,-1,similarity});
                }
            }
            else
            {
                int position_u=-1;
                int position_v=-1;
                int old_similarity=-1;
                for(int i=0;i<(int)neighbor_order[u].size();i++)
                {
                    if(neighbor_order[u][i].second==v)
                    {
                        position_u=i;
                        old_similarity=neighbor_order[u][i].first;
                        break;
                    }
                }
                for(int i=0;i<(int)neighbor_order[v].size();i++)
                {
                    if(neighbor_order[v][i].second==u)
                    {
                        position_v=i;
                        break;
                    }
                }

                int bucket=old_similarity/bucket_width;
                if(bucket==DELTA) bucket=DELTA-1;
                interval_num[u][bucket]--;
                interval_num[v][bucket]--;
                if(bucket>0)
                {
                    affected_vertex_set.insert(u);
                    affected_vertex_set.insert(v);
                    changed_edges.push_back({u,v,old_similarity,-1});
                }
                neighbor_order[u].erase(neighbor_order[u].begin()+position_u);
                neighbor_order[v].erase(neighbor_order[v].begin()+position_v);
            }

            for(auto vertex:dirty_neighbor_order)
            {
                sort(neighbor_order[vertex].begin(),neighbor_order[vertex].end(),cmp_order_space_index);
            }

            affected_vertices.assign(affected_vertex_set.begin(),affected_vertex_set.end());

            // Aggregate CT deltas first. The same vertex may gain and lose
            // different similarity edges at one epsilon, so sequentially
            // changing its core state would expose an invalid intermediate PPT.
            unordered_map<int,vector<int>> vertex_ct_delta;
            for(auto &edge:changed_edges)
            {
                int old_bucket=0;
                int new_bucket=0;
                if(edge.old_similarity>=0)
                {
                    old_bucket=min(edge.old_similarity/bucket_width,DELTA-1);
                }
                if(edge.new_similarity>=0)
                {
                    new_bucket=min(edge.new_similarity/bucket_width,DELTA-1);
                }

                int begin_epsilon=min(old_bucket,new_bucket)+1;
                int end_epsilon=max(old_bucket,new_bucket);
                int change=new_bucket>old_bucket?1:-1;
                int edge_vertices[2]={edge.u,edge.v};
                for(int k=0;k<2;k++)
                {
                    int vertex=edge_vertices[k];
                    if(vertex_ct_delta.find(vertex)==vertex_ct_delta.end())
                    {
                        vertex_ct_delta[vertex]=vector<int>(DELTA,0);
                    }
                    for(int epsilon=begin_epsilon;epsilon<=end_epsilon;epsilon++)
                    {
                        vertex_ct_delta[vertex][epsilon]+=change;
                    }
                }
            }

            for(auto &vertex_delta:vertex_ct_delta)
            {
                int vertex=vertex_delta.first;
                int new_ct=1;
                for(int epsilon=DELTA-1;epsilon>0;epsilon--)
                {
                    new_ct+=interval_num[vertex][epsilon];
                    int delta=vertex_delta.second[epsilon];
                    if(delta==0) continue;

                    int old_ct=new_ct-delta;
                    int begin_mu=min(old_ct,new_ct)+1;
                    int end_mu=max(old_ct,new_ct);
                    affected_mu_intervals[epsilon].push_back(
                        make_pair(begin_mu,end_mu));
                }
            }

            // An active-edge change may alter connectivity even when endpoint
            // CT deltas cancel. Cover every mu for which both endpoints were
            // core in either the old or the final state.
            for(auto &edge:changed_edges)
            {
                int old_bucket=0;
                int new_bucket=0;
                if(edge.old_similarity>=0)
                {
                    old_bucket=min(edge.old_similarity/bucket_width,DELTA-1);
                }
                if(edge.new_similarity>=0)
                {
                    new_bucket=min(edge.new_similarity/bucket_width,DELTA-1);
                }

                int begin_epsilon=min(old_bucket,new_bucket)+1;
                int end_epsilon=max(old_bucket,new_bucket);
                for(int epsilon=begin_epsilon;epsilon<=end_epsilon;epsilon++)
                {
                    int new_ct_u=1;
                    int new_ct_v=1;
                    for(int bucket=epsilon;bucket<DELTA;bucket++)
                    {
                        new_ct_u+=interval_num[edge.u][bucket];
                        new_ct_v+=interval_num[edge.v][bucket];
                    }
                    int delta_u=vertex_ct_delta[edge.u][epsilon];
                    int delta_v=vertex_ct_delta[edge.v][epsilon];
                    int old_ct_u=new_ct_u-delta_u;
                    int old_ct_v=new_ct_v-delta_v;
                    int max_mu=max(min(old_ct_u,old_ct_v),min(new_ct_u,new_ct_v));
                    if(max_mu<2) continue;

                    affected_mu_intervals[epsilon].push_back(
                        make_pair(2,max_mu));
                }
            }

            // Merge only overlapping intervals. Disjoint core and connectivity
            // ranges remain separate so unrelated mu values are not rebuilt.
            for(int epsilon=1;epsilon<DELTA;epsilon++)
            {
                vector<pair<int,int>> &intervals=affected_mu_intervals[epsilon];
                if(intervals.empty()) continue;
                sort(intervals.begin(),intervals.end());
                int write_position=0;
                for(int i=1;i<(int)intervals.size();i++)
                {
                    if(intervals[i].first<=intervals[write_position].second+1)
                    {
                        intervals[write_position].second=
                            max(intervals[write_position].second,intervals[i].second);
                    }
                    else
                    {
                        write_position++;
                        intervals[write_position]=intervals[i];
                    }
                }
                intervals.resize(write_position+1);
            }
        }

        void restore_ppt_membership_for_update()
        {
            free_update_cs_ids.clear();
            point2ppt.clear();
            vertex2ppt.clear();
            vertex2ppt.resize(n);
            vertex2cs.clear();
            vertex2cs.resize(n);
            cs2ppt.clear();
            cs2ppt.resize(cluster_slices.size(),-1);
            delete[] max_mu_under_epsilon;
            max_mu_under_epsilon=new int[DELTA];
            for(int epsilon=0;epsilon<DELTA;epsilon++)
            {
                max_mu_under_epsilon[epsilon]=1;
            }

            // The persistent file stores cluster slices instead of repeating
            // PPT.vertices. Recover both directions and the two-dimensional
            // navigation information before the first update.
            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                PPT &ppt=ppts[ppt_id];
                ppt.id=ppt_id;
                ppt.vertices.clear();
                point2ppt[make_pair(ppt.epsilon,ppt.mu)]=ppt_id;
                max_mu_under_epsilon[ppt.epsilon]=
                    max(max_mu_under_epsilon[ppt.epsilon],ppt.mu);

                for(auto cs_id:ppt.componet_core_ids)
                {
                    for(auto vertex:cluster_slices[cs_id])
                    {
                        ppt.vertices.push_back(vertex);
                        vertex2ppt[vertex].push_back(ppt_id);
                        vertex2cs[vertex].push_back(cs_id);
                    }
                    cs2ppt[cs_id]=ppt_id;
                }
            }

            ppt_matrix.clear();
            ppt_matrix.resize(DELTA);
            for(int epsilon=0;epsilon<DELTA;epsilon++)
            {
                ppt_matrix[epsilon].assign(
                    max_mu_under_epsilon[epsilon]-1,-1);
            }
            for(auto &ppt:ppts)
            {
                ppt_matrix[ppt.epsilon][ppt.mu-2]=ppt.id;
            }
        }

        void rebuild_ppt_maps_after_membership_update()
        {
            vector<PPT> new_ppts;
            new_ppts.reserve(ppts.size());
            for(auto &ppt:ppts)
            {
                if(ppt.vertices.empty())
                {
                    // Reuse slices owned by a disappeared PPT in later updates.
                    for(auto cs_id:ppt.componet_core_ids)
                    {
                        if(cs_id<0||cs_id>=(int)cluster_slices.size()) continue;
                        cluster_slices[cs_id].clear();
                        non_cores[cs_id].clear();
                        free_update_cs_ids.push_back(cs_id);
                    }
                    continue;
                }
                ppt.id=new_ppts.size();
                new_ppts.push_back(move(ppt));
            }
            ppts=move(new_ppts);

            point2ppt.clear();
            vertex2ppt.clear();
            vertex2cs.clear();
            cs2ppt.clear();
            vertex2ppt.resize(n);
            vertex2cs.resize(n);
            cs2ppt.resize(cluster_slices.size(),-1);
            delete[] max_mu_under_epsilon;
            max_mu_under_epsilon=new int[DELTA];
            for(int i=0;i<DELTA;i++) max_mu_under_epsilon[i]=1;

            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                PPT &ppt=ppts[ppt_id];
                ppt.id=ppt_id;
                point2ppt[make_pair(ppt.epsilon,ppt.mu)]=ppt_id;
                max_mu_under_epsilon[ppt.epsilon]=
                    max(max_mu_under_epsilon[ppt.epsilon],ppt.mu);
                for(auto vertex:ppt.vertices)
                {
                    vertex2ppt[vertex].push_back(ppt_id);
                    vertex2cs[vertex].push_back(-1);
                }
                for(auto cs_id:ppt.componet_core_ids)
                {
                    if(cs_id<0||cs_id>=(int)cluster_slices.size()) continue;
                    cs2ppt[cs_id]=ppt_id;
                    for(auto vertex:cluster_slices[cs_id])
                    {
                        for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                        {
                            if(vertex2ppt[vertex][i]==ppt_id)
                            {
                                vertex2cs[vertex][i]=cs_id;
                                break;
                            }
                        }
                    }
                }
            }

            ppt_matrix.clear();
            ppt_matrix.resize(DELTA);
            for(int i=0;i<DELTA;i++)
            {
                ppt_matrix[i].assign(max_mu_under_epsilon[i]-1,-1);
            }
            for(auto &ppt:ppts)
            {
                ppt_matrix[ppt.epsilon][ppt.mu-2]=ppt.id;
            }
        }

        bool update_affected_ppt_membership(
            const vector<int> &affected_vertices,
            vector<pair<int,int>> &membership_changed_points)
        {
            bool membership_changed=false;
            vector<int> touched_ppt_ids;
            vector<char> touched_epsilon(DELTA,0);
            membership_changed_points.clear();
            for(auto vertex:affected_vertices)
            {
                vector<pair<int,int>> old_points;
                vector<pair<int,int>> new_points;
                unordered_map<pair<int,int>,int,PairHash> old_point2cs;
                for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                {
                    int ppt_id=vertex2ppt[vertex][i];
                    pair<int,int> point=make_pair(
                        ppts[ppt_id].epsilon,ppts[ppt_id].mu);
                    int cs_id=-1;
                    if(i<(int)vertex2cs[vertex].size()) cs_id=vertex2cs[vertex][i];
                    old_points.push_back(make_pair(ppts[ppt_id].epsilon,ppts[ppt_id].mu));
                    old_point2cs[point]=cs_id;
                }

                int new_ct=1;
                for(int bucket=DELTA-1;bucket>0;--bucket)
                {
                    new_ct+=interval_num[vertex][bucket];
                    if(interval_num[vertex][bucket]!=0)
                    {
                        new_points.push_back(make_pair(bucket,new_ct));
                    }
                }

                sort(old_points.begin(),old_points.end());
                sort(new_points.begin(),new_points.end());
                if(old_points==new_points) continue;
                membership_changed=true;

                int old_position=0;
                int new_position=0;
                while(old_position<(int)old_points.size()||
                      new_position<(int)new_points.size())
                {
                    if(new_position==(int)new_points.size()||
                       (old_position<(int)old_points.size()&&
                        old_points[old_position]<new_points[new_position]))
                    {
                        membership_changed_points.push_back(old_points[old_position++]);
                    }
                    else if(old_position==(int)old_points.size()||
                            new_points[new_position]<old_points[old_position])
                    {
                        membership_changed_points.push_back(new_points[new_position++]);
                    }
                    else
                    {
                        old_position++;
                        new_position++;
                    }
                }

                // Remove the vertex from every old PPT entry first. Its complete
                // new staircase is then derived from the final interval counts.
                for(auto ppt_id:vertex2ppt[vertex])
                {
                    vector<int> &vertices=ppts[ppt_id].vertices;
                    auto position=find(vertices.begin(),vertices.end(),vertex);
                    if(position!=vertices.end()) vertices.erase(position);
                    touched_ppt_ids.push_back(ppt_id);
                    touched_epsilon[ppts[ppt_id].epsilon]=1;
                }

                vertex2ppt[vertex].clear();
                vertex2cs[vertex].clear();
                for(auto point:new_points)
                {
                    auto ppt_it=point2ppt.find(point);
                    int ppt_id;
                    if(ppt_it==point2ppt.end())
                    {
                        PPT ppt;
                        ppt.id=ppts.size();
                        ppt.epsilon=point.first;
                        ppt.mu=point.second;
                        ppt.vertices.push_back(vertex);
                        ppts.push_back(move(ppt));
                        ppt_id=ppts.size()-1;
                        point2ppt[point]=ppt_id;

                        if((int)ppt_matrix[point.first].size()<point.second-1)
                        {
                            ppt_matrix[point.first].resize(point.second-1,-1);
                        }
                        ppt_matrix[point.first][point.second-2]=ppt_id;
                    }
                    else
                    {
                        ppt_id=ppt_it->second;
                        ppts[ppt_id].vertices.push_back(vertex);
                    }
                    vertex2ppt[vertex].push_back(ppt_id);
                    auto old_cs_it=old_point2cs.find(point);
                    vertex2cs[vertex].push_back(
                        old_cs_it==old_point2cs.end()?-1:old_cs_it->second);
                    touched_ppt_ids.push_back(ppt_id);
                    touched_epsilon[point.first]=1;
                }
            }

            if(membership_changed)
            {
                // Keep PPT IDs stable. Only points touched by changed vertices
                // update the point map and the two-dimensional navigation table.
                sort(touched_ppt_ids.begin(),touched_ppt_ids.end());
                touched_ppt_ids.erase(
                    unique(touched_ppt_ids.begin(),touched_ppt_ids.end()),
                    touched_ppt_ids.end());
                for(auto ppt_id:touched_ppt_ids)
                {
                    PPT &ppt=ppts[ppt_id];
                    pair<int,int> point=make_pair(ppt.epsilon,ppt.mu);
                    if(!ppt.vertices.empty())
                    {
                        point2ppt[point]=ppt_id;
                        ppt_matrix[ppt.epsilon][ppt.mu-2]=ppt_id;
                        continue;
                    }

                    auto point_it=point2ppt.find(point);
                    if(point_it!=point2ppt.end()&&point_it->second==ppt_id)
                    {
                        point2ppt.erase(point_it);
                    }
                    if(ppt.mu-2<(int)ppt_matrix[ppt.epsilon].size()&&
                       ppt_matrix[ppt.epsilon][ppt.mu-2]==ppt_id)
                    {
                        ppt_matrix[ppt.epsilon][ppt.mu-2]=-1;
                    }
                }

                for(int epsilon=1;epsilon<DELTA;epsilon++)
                {
                    if(!touched_epsilon[epsilon]) continue;
                    int max_mu=ppt_matrix[epsilon].size()+1;
                    while(max_mu>=2&&ppt_matrix[epsilon][max_mu-2]==-1)
                    {
                        max_mu--;
                    }
                    max_mu_under_epsilon[epsilon]=max_mu;
                }
            }
            return membership_changed;
        }

        void prepare_affected_component_region(int max_affected_epsilon)
        {
            vector<vector<int>> old_cluster_slices=move(cluster_slices);
            vector<vector<int>> old_non_cores=move(non_cores);
            vector<int> old_cs2new_cs(old_cluster_slices.size(),-1);

            cluster_slices.clear();
            non_cores.clear();
            cs2ppt.clear();

            // Entries above the affected epsilon closure keep their slices.
            // Entries inside the closure discard their component information
            // and will be reconstructed in the original bottom-up order.
            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                PPT &ppt=ppts[ppt_id];
                if(ppt.epsilon<=max_affected_epsilon)
                {
                    ppt.componet_core_ids.clear();
                    ppt.edges.clear();
                    ppt.connectivitys.clear();
                    ppt.cs2edges.clear();
                    continue;
                }

                vector<int> new_component_ids;
                for(auto old_cs_id:ppt.componet_core_ids)
                {
                    int new_cs_id=cluster_slices.size();
                    old_cs2new_cs[old_cs_id]=new_cs_id;
                    new_component_ids.push_back(new_cs_id);
                    cluster_slices.push_back(move(old_cluster_slices[old_cs_id]));
                    non_cores.push_back(move(old_non_cores[old_cs_id]));
                    cs2ppt.push_back(ppt_id);
                }
                ppt.componet_core_ids=move(new_component_ids);

                unordered_map<int,int> new_cs2edges;
                for(auto &item:ppt.cs2edges)
                {
                    int old_cs_id=item.first;
                    if(old_cs_id>=0&&old_cs_id<(int)old_cs2new_cs.size()&&
                       old_cs2new_cs[old_cs_id]!=-1)
                    {
                        new_cs2edges[old_cs2new_cs[old_cs_id]]=item.second;
                    }
                }
                ppt.cs2edges=move(new_cs2edges);
            }

            vertex2cs.clear();
            vertex2cs.resize(n);
            for(int vertex=0;vertex<n;vertex++)
            {
                vertex2cs[vertex].resize(vertex2ppt[vertex].size(),-1);
            }
            for(int cs_id=0;cs_id<(int)cluster_slices.size();cs_id++)
            {
                int ppt_id=cs2ppt[cs_id];
                for(auto vertex:cluster_slices[cs_id])
                {
                    for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                    {
                        if(vertex2ppt[vertex][i]==ppt_id)
                        {
                            vertex2cs[vertex][i]=cs_id;
                            break;
                        }
                    }
                }
            }
        }

        void rebuild_affected_ppt_entries(
            const vector<vector<pair<int,int>>> &affected_mu_intervals,
            const vector<pair<int,int>> &membership_changed_points)
        {
            const auto select_begin=std::chrono::steady_clock().now();
            // Epsilon and mu are both selected precisely. The intervals were
            // merged when changes were collected, but disjoint ranges and their
            // lower bounds must remain intact.
            vector<int> affected_ppt_ids;
            affected_ppt_ids.reserve(ppts.size());
            unordered_set<pair<int,int>,PairHash> membership_changed_set(
                membership_changed_points.begin(),membership_changed_points.end());
            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                PPT &ppt=ppts[ppt_id];
                bool affected=membership_changed_set.find(
                    make_pair(ppt.epsilon,ppt.mu))!=membership_changed_set.end();
                for(auto interval:affected_mu_intervals[ppt.epsilon])
                {
                    if(affected) break;
                    if(ppt.mu>=interval.first&&ppt.mu<=interval.second)
                    {
                        affected=true;
                        break;
                    }
                }
                if(affected) affected_ppt_ids.push_back(ppt_id);
            }
            sort(affected_ppt_ids.begin(),affected_ppt_ids.end(),[this](int a,int b)
            {
                if(ppts[a].epsilon!=ppts[b].epsilon)
                {
                    return ppts[a].epsilon>ppts[b].epsilon;
                }
                return ppts[a].mu>ppts[b].mu;
            });

            vector<char> affected_epsilon(DELTA,0);
            vector<char> affected_ppt(ppts.size(),0);
            for(auto ppt_id:affected_ppt_ids)
            {
                affected_ppt[ppt_id]=1;
                affected_epsilon[ppts[ppt_id].epsilon]=1;
            }
            const auto select_end=std::chrono::steady_clock().now();
            update_profile_select_time+=tms(select_end-select_begin).count();
            update_profile_ppt_num=affected_ppt_ids.size();
            for(int epsilon=1;epsilon<DELTA;epsilon++)
            {
                if(affected_epsilon[epsilon]) update_profile_epsilon_num++;
            }

            // Discard only the output owned by affected entries. The update
            // constructor below follows the original bottom-up code and reads
            // every other entry as finalized dominate-space state.
            const auto cleanup_begin=std::chrono::steady_clock().now();
            for(auto ppt_id:affected_ppt_ids)
            {
                PPT &ppt=ppts[ppt_id];
                for(auto cs_id:ppt.componet_core_ids)
                {
                    if(cs_id<0||cs_id>=(int)cluster_slices.size()) continue;
                    cluster_slices[cs_id].clear();
                    non_cores[cs_id].clear();
                    cs2ppt[cs_id]=-1;
                    free_update_cs_ids.push_back(cs_id);
                }
                for(auto vertex:ppt.vertices)
                {
                    for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                    {
                        if(vertex2ppt[vertex][i]!=ppt_id) continue;
                        vertex2cs[vertex][i]=-1;
                        break;
                    }
                }
                ppt.componet_core_ids.clear();
                ppt.edges.clear();
                ppt.connectivitys.clear();
                ppt.cs2edges.clear();
            }
            const auto cleanup_end=std::chrono::steady_clock().now();
            update_profile_cleanup_time+=tms(cleanup_end-cleanup_begin).count();

            update_componets(affected_epsilon,affected_ppt);
            return;

            vector<vector<int>> ppts_by_epsilon(DELTA);
            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                ppts_by_epsilon[ppts[ppt_id].epsilon].push_back(ppt_id);
            }
            vector<int> core_mark(n,0);
            vector<int> non_core_mark(n,0);
            vector<int> new_core_mark(n,0);
            int core_stamp=0;
            int non_core_stamp=0;
            int new_core_stamp=0;
            int shared_epsilon=-1;
            int previous_mu=INT_MAX;
            vector<int> core_vertices;
            UnionFind core_uf(n);
            UnionFind connectivity_uf(n);

            for(auto ppt_id:affected_ppt_ids)
            {
                PPT &current_ppt=ppts[ppt_id];
                int range=current_ppt.epsilon*(PRECISION/DELTA);
                bool start_new_epsilon=shared_epsilon!=current_ppt.epsilon;

                // Step 1: retire only this entry's old slices. Other entries
                // remain available as finalized dominate-space components.
                for(auto cs_id:current_ppt.componet_core_ids)
                {
                    if(cs_id<0||cs_id>=(int)cluster_slices.size()) continue;
                    cluster_slices[cs_id].clear();
                    non_cores[cs_id].clear();
                    cs2ppt[cs_id]=-1;
                    free_update_cs_ids.push_back(cs_id);
                }
                for(auto vertex:current_ppt.vertices)
                {
                    for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                    {
                        if(vertex2ppt[vertex][i]==ppt_id)
                        {
                            vertex2cs[vertex][i]=-1;
                            break;
                        }
                    }
                }
                current_ppt.componet_core_ids.clear();
                current_ppt.edges.clear();
                current_ppt.connectivitys.clear();
                current_ppt.cs2edges.clear();

                // Step 2: fixed epsilon is processed from large mu to small mu.
                // The core set is monotone, so adjacent affected PPTs share one
                // union-find and only newly admitted vertices scan their edges.
                vector<int> new_core_vertices;
                vector<int> newly_dominating_ppt_ids;
                if(start_new_epsilon)
                {
                    shared_epsilon=current_ppt.epsilon;
                    previous_mu=INT_MAX;
                    core_stamp++;
                    non_core_stamp++;
                    core_vertices.clear();

                    for(int epsilon=current_ppt.epsilon;epsilon<DELTA;epsilon++)
                    {
                        for(auto other_id:ppts_by_epsilon[epsilon])
                        {
                            if(ppts[other_id].mu<current_ppt.mu) break;
                            newly_dominating_ppt_ids.push_back(other_id);
                            for(auto vertex:ppts[other_id].vertices)
                            {
                                if(core_mark[vertex]==core_stamp) continue;
                                core_mark[vertex]=core_stamp;
                                core_vertices.push_back(vertex);
                            }
                        }
                    }
                    core_uf.Init(core_vertices);
                    for(auto vertex:core_vertices)
                    {
                        for(auto &neighbor:neighbor_order[vertex])
                        {
                            if(neighbor.first<range) break;
                            if(core_mark[neighbor.second]==core_stamp&&vertex<neighbor.second)
                            {
                                core_uf.Unite(vertex,neighbor.second);
                            }
                        }
                    }
                }
                else
                {
                    for(int epsilon=current_ppt.epsilon;epsilon<DELTA;epsilon++)
                    {
                        for(auto other_id:ppts_by_epsilon[epsilon])
                        {
                            int other_mu=ppts[other_id].mu;
                            if(other_mu>=previous_mu) continue;
                            if(other_mu<current_ppt.mu) break;
                            newly_dominating_ppt_ids.push_back(other_id);
                            for(auto vertex:ppts[other_id].vertices)
                            {
                                if(core_mark[vertex]==core_stamp) continue;
                                core_mark[vertex]=core_stamp;
                                core_vertices.push_back(vertex);
                                new_core_vertices.push_back(vertex);
                            }
                        }
                    }
                    core_uf.Init(new_core_vertices);
                    for(auto vertex:new_core_vertices)
                    {
                        for(auto &neighbor:neighbor_order[vertex])
                        {
                            if(neighbor.first<range) break;
                            if(core_mark[neighbor.second]==core_stamp)
                            {
                                core_uf.Unite(vertex,neighbor.second);
                            }
                        }
                    }
                }
                previous_mu=current_ppt.mu;
                new_core_stamp++;
                for(auto vertex:new_core_vertices)
                {
                    new_core_mark[vertex]=new_core_stamp;
                }

                // Step 3: vertices owned by this PPT are partitioned according
                // to their components in the complete dominate core graph.
                unordered_map<int,vector<int>> root2vertices;
                root2vertices.reserve(current_ppt.vertices.size());
                for(auto vertex:current_ppt.vertices)
                {
                    root2vertices[core_uf.Find(vertex)].push_back(vertex);
                }
                for(auto &component:root2vertices)
                {
                    int cs_id;
                    if(free_update_cs_ids.empty())
                    {
                        cs_id=cluster_slices.size();
                        cluster_slices.push_back(vector<int>());
                        non_cores.push_back(vector<int>());
                        cs2ppt.push_back(-1);
                    }
                    else
                    {
                        cs_id=free_update_cs_ids.back();
                        free_update_cs_ids.pop_back();
                    }
                    cluster_slices[cs_id]=move(component.second);
                    non_cores[cs_id].clear();
                    cs2ppt[cs_id]=ppt_id;
                    current_ppt.componet_core_ids.push_back(cs_id);

                    for(auto vertex:cluster_slices[cs_id])
                    {
                        for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                        {
                            if(vertex2ppt[vertex][i]==ppt_id)
                            {
                                vertex2cs[vertex][i]=cs_id;
                                break;
                            }
                        }
                    }
                }

                // Step 4: pre_non_cores is also monotone while mu decreases.
                // Inject only newly dominating PPTs; non-cores generated by the
                // current entry remain marked for every following lower mu.
                for(auto other_id:newly_dominating_ppt_ids)
                {
                    if(other_id==ppt_id) continue;
                    for(auto cs_id:ppts[other_id].componet_core_ids)
                    {
                        for(auto vertex:non_cores[cs_id])
                        {
                            non_core_mark[vertex]=non_core_stamp;
                        }
                    }
                }
                for(auto cs_id:current_ppt.componet_core_ids)
                {
                    for(auto vertex:cluster_slices[cs_id])
                    {
                        for(auto &neighbor:neighbor_order[vertex])
                        {
                            if(neighbor.first<range) break;
                            int other=neighbor.second;
                            if(core_mark[other]==core_stamp||
                               non_core_mark[other]==non_core_stamp)
                            {
                                continue;
                            }
                            non_core_mark[other]=non_core_stamp;
                            non_cores[cs_id].push_back(other);
                        }
                    }
                }

                // Step 5: UF2 records connectivity already represented by the
                // index. It is initialized once per epsilon and then advances
                // with UF as mu decreases. Newly dominating unchanged PPTs are
                // injected directly; the current rebuilt PPT contributes its
                // new CS before missing cross-CS edges are selected.
                if(start_new_epsilon)
                {
                    connectivity_uf.Init(core_vertices);
                }
                else
                {
                    connectivity_uf.Init(new_core_vertices);
                }
                for(auto other_id:newly_dominating_ppt_ids)
                {
                    for(auto cs_id:ppts[other_id].componet_core_ids)
                    {
                        if(cluster_slices[cs_id].empty()) continue;
                        int first_vertex=cluster_slices[cs_id][0];
                        for(int i=1;i<(int)cluster_slices[cs_id].size();i++)
                        {
                            connectivity_uf.Unite(first_vertex,cluster_slices[cs_id][i]);
                        }
                    }
                    for(int edge_id=0;edge_id<(int)ppts[other_id].edges.size();edge_id++)
                    {
                        if(ppts[other_id].connectivitys[edge_id]<range) continue;
                        auto &edge=ppts[other_id].edges[edge_id];
                        connectivity_uf.Unite(edge.first,edge.second);
                    }
                }

                vector<int> &candidate_vertices=
                    start_new_epsilon?core_vertices:new_core_vertices;
                for(auto vertex:candidate_vertices)
                {
                    for(auto &neighbor:neighbor_order[vertex])
                    {
                        if(neighbor.first<range) break;
                        int other=neighbor.second;
                        if(core_mark[other]!=core_stamp) continue;
                        if(start_new_epsilon&&vertex>=other) continue;
                        if(!start_new_epsilon&&
                           new_core_mark[other]==new_core_stamp&&vertex>=other) continue;
                        if(connectivity_uf.Find(vertex)==connectivity_uf.Find(other)) continue;

                        connectivity_uf.Unite(vertex,other);
                        current_ppt.edges.push_back(make_pair(vertex,other));
                        current_ppt.connectivitys.push_back(neighbor.first);
                    }
                }
            }
        }

        void update_componets(const vector<char> &affected_epsilon,
                              const vector<char> &affected_ppt)
        {
            vector<int> min_affected_mu(DELTA,INF);
            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                if(!affected_ppt[ppt_id]) continue;
                min_affected_mu[ppts[ppt_id].epsilon]=min(
                    min_affected_mu[ppts[ppt_id].epsilon],ppts[ppt_id].mu);
            }

            vector<int> pre_subgraph;
            vector<int> pre_non_cores;
            vector<bool> vertex_bit_set(n,false);
            vector<bool> non_core_bit_set(n,false);
            vector<int> vertex2currentPPT(n,-1);
            vector<int> vertex2currentcs(n,-1);
            UnionFind UF(n);
            UnionFind UF2(n);

            pre_subgraph.reserve(n);
            pre_non_cores.reserve(n);

            // Keep the original epsilon/mu bottom-up construction order, but
            // enter only epsilon values containing an affected PPT entry.
            for(int temp_eps=DELTA-1;temp_eps>0;--temp_eps)
            {
                if(!affected_epsilon[temp_eps]||ppt_matrix[temp_eps].empty()) continue;

                const auto epsilon_state_begin=std::chrono::steady_clock().now();
                int temp_range=temp_eps*(PRECISION/DELTA);
                UF.Init(pre_subgraph);
                UF2.Init(pre_subgraph);
                for(auto vertex:pre_subgraph)
                {
                    vertex_bit_set[vertex]=false;
                    vertex2currentPPT[vertex]=-1;
                    vertex2currentcs[vertex]=-1;
                }
                for(auto vertex:pre_non_cores)
                {
                    non_core_bit_set[vertex]=false;
                }
                pre_subgraph.clear();
                pre_non_cores.clear();
                const auto epsilon_state_end=std::chrono::steady_clock().now();
                update_profile_state_time+=
                    tms(epsilon_state_end-epsilon_state_begin).count();

                // Higher mu values establish the bottom-up state. Levels below
                // the lowest affected PPT do not contribute to rebuilt entries.
                // PPT entries and ppt_matrix only represent mu >= 2.
                int min_update_mu=max(2,min_affected_mu[temp_eps]);
                for(int temp_mu=max_mu_under_epsilon[temp_eps];
                    temp_mu>=min_update_mu;--temp_mu)
                {
                    update_profile_mu_num++;
                    const auto state_begin=std::chrono::steady_clock().now();
                    vector<int> subgraph;
                    vector<int> dominate_space;
                    int current_ppt_id=ppt_matrix[temp_eps][temp_mu-2];
                    if(current_ppt_id!=-1) dominate_space.push_back(current_ppt_id);

                    for(int larger_eps=temp_eps+1;larger_eps<DELTA;larger_eps++)
                    {
                        if(max_mu_under_epsilon[larger_eps]<temp_mu) continue;
                        int dominate_id=ppt_matrix[larger_eps][temp_mu-2];
                        if(dominate_id!=-1) dominate_space.push_back(dominate_id);
                    }

                    if(temp_mu==max_mu_under_epsilon[temp_eps])
                    {
                        for(int larger_eps=temp_eps+1;larger_eps<DELTA;larger_eps++)
                        {
                            for(int larger_mu=temp_mu+1;
                                larger_mu<=max_mu_under_epsilon[larger_eps];larger_mu++)
                            {
                                int dominate_id=ppt_matrix[larger_eps][larger_mu-2];
                                if(dominate_id!=-1) dominate_space.push_back(dominate_id);
                            }
                        }
                    }

                    // This is the original incremental subgraph construction.
                    for(auto dominate_id:dominate_space)
                    {
                        for(auto vertex:ppts[dominate_id].vertices)
                        {
                            if(vertex_bit_set[vertex]) continue;
                            vertex_bit_set[vertex]=true;
                            vertex2currentPPT[vertex]=dominate_id;
                            pre_subgraph.push_back(vertex);
                            subgraph.push_back(vertex);
                        }
                        for(auto cs_id:ppts[dominate_id].componet_core_ids)
                        {
                            for(auto vertex:non_cores[cs_id])
                            {
                                if(non_core_bit_set[vertex]) continue;
                                non_core_bit_set[vertex]=true;
                                pre_non_cores.push_back(vertex);
                            }
                        }
                    }
                    update_profile_subgraph_vertex_num+=subgraph.size();
                    const auto state_end=std::chrono::steady_clock().now();
                    update_profile_state_time+=tms(state_end-state_begin).count();

                    const auto uf1_begin=std::chrono::steady_clock().now();
                    for(auto vertex:subgraph)
                    {
                        for(auto &neighbor:neighbor_order[vertex])
                        {
                            if(neighbor.first<temp_range) break;
                            update_profile_uf1_neighbor_num++;
                            if(vertex_bit_set[neighbor.second]&&
                               UF.Find(vertex)!=UF.Find(neighbor.second))
                            {
                                UF.Unite(vertex,neighbor.second);
                            }
                        }
                    }
                    const auto uf1_end=std::chrono::steady_clock().now();
                    update_profile_uf1_time+=tms(uf1_end-uf1_begin).count();

                    // Only an affected real PPT produces new CS and non-cores.
                    const auto component_begin=std::chrono::steady_clock().now();
                    if(current_ppt_id!=-1&&affected_ppt[current_ppt_id])
                    {
                        unordered_map<int,vector<int>> root2vertices;
                        update_profile_cs_vertex_num+=ppts[current_ppt_id].vertices.size();
                        for(auto vertex:ppts[current_ppt_id].vertices)
                        {
                            root2vertices[UF.Find(vertex)].push_back(vertex);
                        }
                        for(auto &component:root2vertices)
                        {
                            int cs_id;
                            if(free_update_cs_ids.empty())
                            {
                                cs_id=cluster_slices.size();
                                cluster_slices.push_back(vector<int>());
                                non_cores.push_back(vector<int>());
                                cs2ppt.push_back(-1);
                            }
                            else
                            {
                                cs_id=free_update_cs_ids.back();
                                free_update_cs_ids.pop_back();
                            }
                            cluster_slices[cs_id]=move(component.second);
                            non_cores[cs_id].clear();
                            cs2ppt[cs_id]=current_ppt_id;
                            ppts[current_ppt_id].componet_core_ids.push_back(cs_id);

                            for(auto vertex:cluster_slices[cs_id])
                            {
                                for(int i=0;i<(int)vertex2ppt[vertex].size();i++)
                                {
                                    if(vertex2ppt[vertex][i]!=current_ppt_id) continue;
                                    vertex2cs[vertex][i]=cs_id;
                                    break;
                                }
                            }
                        }

                        for(auto cs_id:ppts[current_ppt_id].componet_core_ids)
                        {
                            for(auto vertex:cluster_slices[cs_id])
                            {
                                for(auto &neighbor:neighbor_order[vertex])
                                {
                                    if(neighbor.first<temp_range) break;
                                    update_profile_non_core_neighbor_num++;
                                    int other=neighbor.second;
                                    if(vertex_bit_set[other]||non_core_bit_set[other]) continue;
                                    non_core_bit_set[other]=true;
                                    pre_non_cores.push_back(other);
                                    non_cores[cs_id].push_back(other);
                                }
                            }
                        }
                    }
                    const auto component_end=std::chrono::steady_clock().now();
                    update_profile_component_time+=tms(component_end-component_begin).count();

                    // Existing unaffected PPTs and newly rebuilt PPTs both
                    // contribute their stored representation to UF2.
                    const auto uf2_begin=std::chrono::steady_clock().now();
                    for(auto dominate_id:dominate_space)
                    {
                        for(auto cs_id:ppts[dominate_id].componet_core_ids)
                        {
                            if(cluster_slices[cs_id].empty()) continue;
                            update_profile_uf2_cs_vertex_num+=cluster_slices[cs_id].size();
                            int first_vertex=cluster_slices[cs_id][0];
                            vertex2currentcs[first_vertex]=cs_id;
                            for(int i=1;i<(int)cluster_slices[cs_id].size();i++)
                            {
                                vertex2currentcs[cluster_slices[cs_id][i]]=cs_id;
                                UF2.Unite(first_vertex,cluster_slices[cs_id][i]);
                            }
                        }
                        update_profile_uf2_edge_num+=ppts[dominate_id].edges.size();
                        for(int edge_id=0;
                            edge_id<(int)ppts[dominate_id].edges.size();edge_id++)
                        {
                            if(ppts[dominate_id].connectivitys[edge_id]<temp_range) continue;
                            auto &edge=ppts[dominate_id].edges[edge_id];
                            UF2.Unite(edge.first,edge.second);
                        }
                    }
                    const auto uf2_end=std::chrono::steady_clock().now();
                    update_profile_uf2_time+=tms(uf2_end-uf2_begin).count();

                    const auto connectivity_begin=std::chrono::steady_clock().now();
                    for(auto vertex:subgraph)
                    {
                        int owner_ppt_id=vertex2currentPPT[vertex];
                        for(auto &neighbor:neighbor_order[vertex])
                        {
                            if(neighbor.first<temp_range) break;
                            update_profile_connectivity_neighbor_num++;
                            int other=neighbor.second;
                            if(!vertex_bit_set[other]||
                               UF2.Find(vertex)==UF2.Find(other)) continue;

                            // UF2 already contains every retained index edge.
                            // A remaining disconnection therefore needs a new
                            // persistent connectivity edge. Store it under the
                            // original bottom-up owner even when that PPT did
                            // not need its CS/non-core entry reconstructed.
                            UF2.Unite(vertex,other);
                            int other_cs_id=vertex2currentcs[other];
                            auto edge_it=ppts[owner_ppt_id].cs2edges.find(other_cs_id);
                            if(edge_it==ppts[owner_ppt_id].cs2edges.end())
                            {
                                update_profile_new_edge_num++;
                                ppts[owner_ppt_id].cs2edges[other_cs_id]=
                                    ppts[owner_ppt_id].edges.size();
                                ppts[owner_ppt_id].edges.push_back(make_pair(vertex,other));
                                ppts[owner_ppt_id].connectivitys.push_back(neighbor.first);
                            }
                            else
                            {
                                int edge_id=edge_it->second;
                                if(neighbor.first<=ppts[owner_ppt_id].connectivitys[edge_id]) continue;
                                update_profile_replace_edge_num++;
                                ppts[owner_ppt_id].edges[edge_id]=make_pair(vertex,other);
                                ppts[owner_ppt_id].connectivitys[edge_id]=neighbor.first;
                            }
                        }
                    }
                    const auto connectivity_end=std::chrono::steady_clock().now();
                    update_profile_connectivity_time+=
                        tms(connectivity_end-connectivity_begin).count();
                }
            }
        }

        vector<vector<int>> collect_clusters_for_update_check(
            int bucket,int mu,bool include_non_cores)
        {
            vector<int> selected_ppts;
            vector<int> selected_cs;
            vector<int> core_mark(n,0);
            vector<int> core_vertices;
            UnionFind check_uf(n);

            for(int ppt_id=0;ppt_id<(int)ppts.size();ppt_id++)
            {
                if(ppts[ppt_id].epsilon<bucket||ppts[ppt_id].mu<mu) continue;
                selected_ppts.push_back(ppt_id);
                for(auto cs_id:ppts[ppt_id].componet_core_ids)
                {
                    if(cluster_slices[cs_id].empty()) continue;
                    selected_cs.push_back(cs_id);
                    for(auto vertex:cluster_slices[cs_id])
                    {
                        if(core_mark[vertex]) continue;
                        core_mark[vertex]=1;
                        core_vertices.push_back(vertex);
                    }
                }
            }

            check_uf.Init(core_vertices);
            // A core vertex can occur in slices of several selected PPTs. Use
            // vertex IDs directly so this check is independent of whichever
            // slice representative happens to be visited last by the query.
            for(auto cs_id:selected_cs)
            {
                int representative=cluster_slices[cs_id][0];
                for(int i=1;i<(int)cluster_slices[cs_id].size();i++)
                {
                    check_uf.Unite(representative,cluster_slices[cs_id][i]);
                }
            }
            int range=bucket*(PRECISION/DELTA);
            for(auto ppt_id:selected_ppts)
            {
                for(int edge_id=0;edge_id<(int)ppts[ppt_id].edges.size();edge_id++)
                {
                    if(ppts[ppt_id].connectivitys[edge_id]<range) continue;
                    int u=ppts[ppt_id].edges[edge_id].first;
                    int v=ppts[ppt_id].edges[edge_id].second;
                    if(!core_mark[u]||!core_mark[v]) continue;
                    check_uf.Unite(u,v);
                }
            }

            unordered_map<int,vector<int>> root2cluster;
            for(auto cs_id:selected_cs)
            {
                int root=check_uf.Find(cluster_slices[cs_id][0]);
                vector<int> &cluster=root2cluster[root];
                cluster.insert(cluster.end(),cluster_slices[cs_id].begin(),cluster_slices[cs_id].end());
                if(include_non_cores)
                {
                    cluster.insert(cluster.end(),non_cores[cs_id].begin(),non_cores[cs_id].end());
                }
            }

            vector<vector<int>> clusters;
            for(auto &item:root2cluster)
            {
                vector<int> &cluster=item.second;
                sort(cluster.begin(),cluster.end());
                cluster.erase(unique(cluster.begin(),cluster.end()),cluster.end());
                clusters.push_back(move(cluster));
            }
            sort(clusters.begin(),clusters.end());
            return clusters;
        }

        vector<int> collect_non_cores_for_update_check(int bucket,int mu)
        {
            vector<int> result;
            vector<int> core_mark(n,0);
            for(auto &ppt:ppts)
            {
                if(ppt.epsilon<bucket||ppt.mu<mu) continue;
                for(auto cs_id:ppt.componet_core_ids)
                {
                    for(auto vertex:cluster_slices[cs_id])
                    {
                        core_mark[vertex]=1;
                    }
                }
            }
            for(auto &ppt:ppts)
            {
                if(ppt.epsilon<bucket||ppt.mu<mu) continue;
                for(auto cs_id:ppt.componet_core_ids)
                {
                    for(auto vertex:non_cores[cs_id])
                    {
                        // Query processing labels all selected core vertices
                        // before considering non-cores. Stale/redundant border
                        // occurrences therefore have no semantic effect.
                        if(!core_mark[vertex]) result.push_back(vertex);
                    }
                }
            }
            sort(result.begin(),result.end());
            result.erase(unique(result.begin(),result.end()),result.end());
            return result;
        }

        void verify_ppt_update_by_full_reconstruction()
        {
            vector<vector<vector<vector<int>>>> incremental_results(
                DELTA,vector<vector<vector<int>>>(dmax+2));
            vector<vector<vector<int>>> incremental_non_core_sets(
                DELTA,vector<vector<int>>(dmax+2));
            for(int bucket=1;bucket<DELTA;bucket++)
            {
                for(int mu=2;mu<=dmax+1;mu++)
                {
                    incremental_results[bucket][mu]=
                        collect_clusters_for_update_check(bucket,mu,false);
                    incremental_non_core_sets[bucket][mu]=
                        collect_non_cores_for_update_check(bucket,mu);
                }
            }

            // Keep the incremental index because subsequent updates must continue
            // from the state produced by way 1 rather than from the oracle.
            vector<PPT> incremental_ppts=move(ppts);
            vector<vector<int>> incremental_cluster_slices=move(cluster_slices);
            vector<vector<int>> incremental_non_cores=move(non_cores);
            unordered_map<pair<int,int>,int,PairHash> incremental_point2ppt=move(point2ppt);
            vector<vector<int>> incremental_vertex2ppt=move(vertex2ppt);
            vector<vector<int>> incremental_vertex2cs=move(vertex2cs);
            vector<int> incremental_cs2ppt=move(cs2ppt);
            vector<int> incremental_free_cs_ids=move(free_update_cs_ids);
            vector<vector<int>> incremental_ppt_matrix=move(ppt_matrix);
            vector<vector<int>> incremental_table_ppts=move(table_ppts);
            vector<int> incremental_max_mu(DELTA,1);
            for(int i=0;i<DELTA;i++) incremental_max_mu[i]=max_mu_under_epsilon[i];
            int incremental_ppt_max_id=ppt_max_id;

            double old_total_build_time=total_build_time;
            double old_total_build_space=total_build_space;
            double old_core_vertices_space=core_vertices_space;
            double old_non_core_vertices_space=non_core_vertices_space;
            double old_edge_space=edge_space;
            double old_connectivity_space=connectivity_space;

            clear_ppt_structure_for_update();
            construct_ppt_vertex_map();
            construct_componets();
            prepare_query_structure_after_update();

            bool correct=true;
            bool non_core_error=false;
            int wrong_bucket=-1;
            int wrong_mu=-1;
            for(int bucket=1;bucket<DELTA&&correct;bucket++)
            {
                for(int mu=2;mu<=dmax+1;mu++)
                {
                    vector<vector<int>> rebuilt_result=
                        collect_clusters_for_update_check(bucket,mu,false);
                    if(incremental_results[bucket][mu]!=rebuilt_result)
                    {
                        correct=false;
                        wrong_bucket=bucket;
                        wrong_mu=mu;
                        break;
                    }
                    vector<int> rebuilt_non_cores=
                        collect_non_cores_for_update_check(bucket,mu);
                    if(incremental_non_core_sets[bucket][mu]!=rebuilt_non_cores)
                    {
                        correct=false;
                        non_core_error=true;
                        wrong_bucket=bucket;
                        wrong_mu=mu;
                        break;
                    }
                }
            }

            clear_ppt_structure_for_update();
            ppts=move(incremental_ppts);
            cluster_slices=move(incremental_cluster_slices);
            non_cores=move(incremental_non_cores);
            point2ppt=move(incremental_point2ppt);
            vertex2ppt=move(incremental_vertex2ppt);
            vertex2cs=move(incremental_vertex2cs);
            cs2ppt=move(incremental_cs2ppt);
            free_update_cs_ids=move(incremental_free_cs_ids);
            ppt_matrix=move(incremental_ppt_matrix);
            table_ppts=move(incremental_table_ppts);
            max_mu_under_epsilon=new int[DELTA];
            for(int i=0;i<DELTA;i++) max_mu_under_epsilon[i]=incremental_max_mu[i];
            ppt_max_id=incremental_ppt_max_id;

            total_build_time=old_total_build_time;
            total_build_space=old_total_build_space;
            core_vertices_space=old_core_vertices_space;
            non_core_vertices_space=old_non_core_vertices_space;
            edge_space=old_edge_space;
            connectivity_space=old_connectivity_space;

            if(correct)
            {
                cout<<"[UPDATE CHECK] PPT way 1 passed"<<endl;
            }
            else
            {
                cout<<"[UPDATE CHECK ERROR] epsilon bucket "<<wrong_bucket
                    <<" mu "<<wrong_mu
                    <<(non_core_error?" non-core coverage":" core connectivity")<<endl;
            }
        }

        void update(string update_type, int u, int v, int update_way) override
        {
            const auto similarity_begin=std::chrono::steady_clock().now();
            if(neighbor_order.empty())
            {
                fast_construct_neighbor_order();
            }

            if((update_way==1||update_way==3)&&ppts.empty())
            {
                load_index();
            }

            vector<int> affected_vertices;
            vector<vector<pair<int,int>>> affected_mu_intervals;
            update_similarity_graph_for_ppt(
                update_type,u,v,affected_vertices,affected_mu_intervals);
            const auto similarity_end=std::chrono::steady_clock().now();
            tms similarity_time=similarity_end-similarity_begin;

            if(update_way==0)
            {
                // Rebuild the complete PPT layer from the final similarity
                // graph, matching Forest way 0's reconstruction baseline.
                const auto update_begin=std::chrono::steady_clock().now();
                clear_ppt_structure_for_update();
                construct_ppt_vertex_map();
                construct_componets();
                prepare_query_structure_after_update();

                const auto update_end=std::chrono::steady_clock().now();
                tms update_time=update_end-update_begin;
                total_update_time+=update_time.count();
                update_count++;
                cout<<"[UPDATE] update time : "<<update_time.count()<<endl;
                return;
            }
            else if(update_way==1||update_way==3)
            {
                // Recompute affected PPT memberships, then rebuild their safe
                // downward dominance closure. Higher epsilon entries retain
                // their old slices and do not repeat component construction.
                const auto update_begin=std::chrono::steady_clock().now();
                update_profile_select_time=0.0;
                update_profile_cleanup_time=0.0;
                update_profile_state_time=0.0;
                update_profile_uf1_time=0.0;
                update_profile_component_time=0.0;
                update_profile_uf2_time=0.0;
                update_profile_connectivity_time=0.0;
                update_profile_epsilon_num=0;
                update_profile_ppt_num=0;
                update_profile_mu_num=0;
                update_profile_subgraph_vertex_num=0;
                update_profile_uf1_neighbor_num=0;
                update_profile_cs_vertex_num=0;
                update_profile_non_core_neighbor_num=0;
                update_profile_uf2_cs_vertex_num=0;
                update_profile_uf2_edge_num=0;
                update_profile_connectivity_neighbor_num=0;
                update_profile_new_edge_num=0;
                update_profile_replace_edge_num=0;
                update_profile_edge_num_before=0;
                update_profile_edge_num_after=0;
                for(auto &ppt:ppts)
                {
                    update_profile_edge_num_before+=ppt.edges.size();
                }
                vector<pair<int,int>> membership_changed_points;
                const auto membership_begin=std::chrono::steady_clock().now();
                update_affected_ppt_membership(
                    affected_vertices,membership_changed_points);
                const auto membership_end=std::chrono::steady_clock().now();
                tms membership_time=membership_end-membership_begin;
                if(affected_vertices.empty())
                {
                    const auto update_end=std::chrono::steady_clock().now();
                    tms update_time=update_end-update_begin;
                    total_update_time+=update_time.count();
                    update_count++;
                    cout<<"[UPDATE] PPT membership does not change"<<endl;
                    if(update_way==3)
                    {
                        verify_ppt_update_by_full_reconstruction();
                    }
                    cout<<"[UPDATE] update time : "<<update_time.count()<<endl;
                    return;
                }
                table_ppts.clear();
                core_vertices_space=0.0;
                non_core_vertices_space=0.0;
                edge_space=0.0;
                connectivity_space=0.0;
                total_build_space=0.0;
                rebuild_affected_ppt_entries(
                    affected_mu_intervals,membership_changed_points);
                for(auto &ppt:ppts)
                {
                    update_profile_edge_num_after+=ppt.edges.size();
                }
                const auto query_begin=std::chrono::steady_clock().now();
                prepare_query_structure_after_update();
                const auto query_end=std::chrono::steady_clock().now();
                tms query_time=query_end-query_begin;
                const auto update_end=std::chrono::steady_clock().now();

                if(update_way==3)
                {
                    verify_ppt_update_by_full_reconstruction();
                }

                tms update_time=update_end-update_begin;
                total_update_time+=update_time.count();
                update_count++;
                cout<<"[UPDATE] update time : "<<update_time.count()<<endl;
                cout<<"[UPDATE PROFILE] similarity: "<<similarity_time.count()
                    <<" membership/maps: "<<membership_time.count()
                    <<" select affected PPT: "<<update_profile_select_time
                    <<" cleanup entries: "<<update_profile_cleanup_time<<endl;
                cout<<"[UPDATE PROFILE] component state: "<<update_profile_state_time
                    <<" UF1 connectivity: "<<update_profile_uf1_time
                    <<" rebuild CS/non-core: "<<update_profile_component_time
                    <<" inject UF2: "<<update_profile_uf2_time
                    <<" search connectivity: "<<update_profile_connectivity_time<<endl;
                cout<<"[UPDATE PROFILE] query structure: "<<query_time.count()
                    <<" affected epsilon: "<<update_profile_epsilon_num
                    <<" affected PPT: "<<update_profile_ppt_num
                    <<" processed mu: "<<update_profile_mu_num
                    <<" added subgraph vertices: "<<update_profile_subgraph_vertex_num
                    <<endl;
                cout<<"[UPDATE WORK] UF1 neighbor visits: "
                    <<update_profile_uf1_neighbor_num
                    <<" CS vertices: "<<update_profile_cs_vertex_num
                    <<" non-core neighbor visits: "
                    <<update_profile_non_core_neighbor_num<<endl;
                cout<<"[UPDATE WORK] UF2 CS vertices: "
                    <<update_profile_uf2_cs_vertex_num
                    <<" UF2 edges: "<<update_profile_uf2_edge_num
                    <<" connectivity neighbor visits: "
                    <<update_profile_connectivity_neighbor_num<<endl;
                cout<<"[UPDATE WORK] new edges: "<<update_profile_new_edge_num
                    <<" replaced edges: "<<update_profile_replace_edge_num
                    <<" stored edges before: "<<update_profile_edge_num_before
                    <<" after: "<<update_profile_edge_num_after<<endl;
                return;
            }
        }
        void print_update_time() override
        {
            if(update_count==0) return;
            cout<<"[PRINT TIME] average update time: "
                <<total_update_time/double(update_count)<<endl;
        }    
};
