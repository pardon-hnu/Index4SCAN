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
        
        int * max_mu_under_epsilon; // when fix a certain epsilon, the max mu that (epsilon,mu) is in space.
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
        void construct_componets()
        {
            cout<<"[CONSTRUCT INDEX] start building componets"<<endl;
            const auto begin=std::chrono::steady_clock().now();
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

            for(int temp_eps=DELTA-1;temp_eps>=0;--temp_eps)
            {
                //if doesn't exist
                if(!ppt_matrix[temp_eps].size()) continue;
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
                            if(vertex_bit_set[nv.second]&&UF.Find(v)!=UF.Find(nv.second))
                            {UF.Unite(v,nv.second);}
                        }
                    }

                    //divide componet for this skyline parameter pair
                    if(current_ppt_id!=-1)
                    {
                        unordered_map<int,vector<int>> root2vertex;
                        for(auto & v:ppts[current_ppt_id].vertices)
                        {root2vertex[UF.Find(v)].push_back(v);}
                        for(auto & item:root2vertex)
                        {
                            ppts[current_ppt_id].componet_core_ids.push_back(cluster_slices.size());
                            cluster_slices.push_back(item.second);
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
                                        ppts[v_ppt_id].edges[edge_id]=make_pair(v,nv.second);
                                        ppts[v_ppt_id].connectivitys[edge_id]=nv.first;
                                    }
                                }
                                else
                                {
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
};