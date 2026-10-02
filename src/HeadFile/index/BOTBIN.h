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
#include <util/HashSet.h>
#include <index/Index.h>
#include <index/GS_Index.h>
#include <set>
#include <memory>
using namespace std;
typedef std::chrono::duration<double> tms;
bool cmp_order_botbin(pair<int,int> x,pair<int,int> y)
{
	return x.first>y.first;
}
bool cmp_bucket(pair<int,int> x,pair<int,int> y)
{
	return x.first>y.first;
}
class BOTBIN : public SCAN_Index
{
    private:
        int k=0;
        float RHO=0.1;
        float FAILURE_PB=0.001;
        int DELTA=10;
        static vector<int> rnk;
        vector<vector<int> > sketch;
        vector<vector<pair<int, int>>> bucket;
        vector<vector<pair<int, int>>> neighbor_order;
        double total_ari=0.0;
        int query_count=0;
        double bucket_space=0.0;
        double neighbor_order_space=0.0;
        double graph_space=0.0;
        double total_bfs_time=0.0;
    public:
        BOTBIN(string dataset,float rho=0.1,float failure_pb=0.001,int delta=10)
        {
            this->RHO=rho;
            this->FAILURE_PB=failure_pb;
            this->DELTA=delta;
            k=(int)(ceil((1/(2*pow(RHO,2)))*(log(2/FAILURE_PB))));
            init_graph(dataset);
        }
        ~BOTBIN()
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
        int locate(int epsilon)
        {
            int bucket_range=PRECISION/DELTA;
            int bucket_id= epsilon/bucket_range;
            return bucket_id;
        }
        void cluster_bobtin(int epsilon,int mu)
        {
            clear_cluster_result();
            cout<<"[CLUSTER] clustering with BOTBIN start !!!"<<endl;
            double query_bfs_time=0.0;
            const auto begin=std::chrono::steady_clock().now();
            int locate_bucket=locate(epsilon);
            cout<<"search from bucket["<<locate_bucket<<"]"<<endl;
            for(auto &core:bucket[locate_bucket])
            {
                if(labels[core.second]!=-1) continue;
                if(core.first<mu-1) break;
                // cout<<"begin cluster from core: "<<core.second<<endl;
                const auto bfs_begin=std::chrono::steady_clock().now();
                queue<int> cluster;
                unordered_set<int> cluster_set;
                cluster.push(core.second);
                cluster_set.insert(core.second);
                // cout<<"cluster path: ";
                while(!cluster.empty())
                {
                    int v=cluster.front();
                    cluster.pop();
                    // cout<<v<<" ";
                    labels[v]=cluster_number;
                    for(auto &nv:neighbor_order[v])
                    {
                        if(labels[nv.second]!=-1) continue;
                        if(nv.first<epsilon) break;
                        if(neighbor_order[nv.second].size()>=mu-1&&neighbor_order[nv.second][mu-2].first>=epsilon) //core
                        {
                            if(cluster_set.find(nv.second)==cluster_set.end())
                            {
                                cluster.push(nv.second);
                                cluster_set.insert(nv.second);
                            }
                        }
                        else //non core
                        {
                            labels[nv.second]=cluster_number;
                        }
                    }
                }
                // cout<<endl;
                cluster_number++;
                const auto bfs_end=std::chrono::steady_clock().now();
                tms bfs_time=bfs_end-bfs_begin;
                // cout<<"[TIME COST] "<<bfs_time.count()<<" s."<<endl;
                // cout<<"cluster time: "<<bfs_time.count()<<" s."<<endl;
                query_bfs_time+=bfs_time.count();
            }
            const auto end=std::chrono::steady_clock().now();
            tms cluster_time=end-begin;
            total_bfs_time+=query_bfs_time;
            cout<<"[QUERY PROFILE] epsilon="<<double(epsilon)/PRECISION<<" mu="<<mu
                <<" query_time_s="<<cluster_time.count()<<" bfs_time_s="<<query_bfs_time
                <<" bfs_ratio="<<(cluster_time.count()>0.0 ? query_bfs_time/cluster_time.count() : 0.0)<<endl;
            cout<<"[CLSUTER] finish clustering!!!"<<endl;
            cout<<"[TIME COST] "<<cluster_time.count()<<" s."<<endl;
            total_cluster_time+=cluster_time.count();
            query_count++;
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
                sort(neighbor_order[u].begin(),neighbor_order[u].end(),cmp_order_botbin);
            }
            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish building Neighbor Order"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            total_build_time+=construct_time.count();
            for(auto &item:neighbor_order)
            {
                total_build_space+=item.size()*8;
            }
            neighbor_order_space=total_build_space;
            graph_space=total_build_space/2;
        }
        void construct_bucket_index()
        {
            cout<<"[CONSTRUCT INDEX] start building bucket index"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            bucket.resize(DELTA);
            for(int v=0;v<n;v++)
            {
                for(int i=0;i<DELTA;i++)
                {
                    int similairty_range=i*(PRECISION/DELTA);
                    int cn=0;
                    for(auto &w:neighbor_order[v])
                    {
                        if(w.first>=similairty_range)
                        {
                            cn++;
                        }
                        else
                        {
                            break;
                        }
                    }
                    if(cn) bucket[i].push_back(make_pair(cn,v));
                }
            }
            for(int i=0;i<DELTA;i++)
            {
                sort(bucket[i].begin(),bucket[i].end(),cmp_bucket);
            }
            const auto end=std::chrono::steady_clock().now();
            tms construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish building bucket index"<<endl;
            cout<<"[TIME COST] "<<construct_time.count()<<" s."<<endl;
            total_build_time+=construct_time.count();
            for(auto &item:bucket)
            {
                total_build_space+=item.size()*8;
            }
            bucket_space=total_build_space-neighbor_order_space;
        }
        void construct() override
        {
            construct_bottom_k_sketch();
            construct_neighbor_order();
            construct_bucket_index();
        }
        void obtain_space_ratio()
        {
            for(auto &item:neighbor_order)
            {
                neighbor_order_space+=item.size()*8;
            }
            graph_space=neighbor_order_space/2;
            for(auto &item:bucket)
            {
                bucket_space+=item.size()*8;
            }
            cout<<"[PRINT SPACE INFORMATION] neighbor order space cost: "<<(neighbor_order_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] bucket index space cost: "<<(bucket_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] graph space cost: "<<(graph_space/1024)/1024<<".MB"<<endl;
        }
        void query(float epsilon, int mu) override
        {
            int eps=epsilon*PRECISION;
            cluster_bobtin(eps,mu);
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
        void print_build_cost() override
        {
            cout<<"[PRINT TIME] total building time: "<<total_build_time<<endl;
            cout<<"[PRINT SPAEC]: total space cost: "<<(total_build_space/1024)/1024<<" MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] neighbor order space cost: "<<(neighbor_order_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] bucket index space cost: "<<(bucket_space/1024)/1024<<".MB"<<endl;
            cout<<"[PRINT SPACE INFORMATION] graph space cost: "<<(graph_space/1024)/1024<<".MB"<<endl;
        }
        void print_cluster_time() override
        {
            cout<<"[PRINT TIME] total clustering time: "<<total_cluster_time<<endl;
            cout<<"[PRINT TIME] total bfs time: "<<total_bfs_time<<endl;
            cout<<"[PRINT TIME] average bfs time: "<<(query_count>0 ? total_bfs_time/query_count : 0.0)<<endl;
            cout<<"[PRINT TIME] bfs ratio: "<<(total_cluster_time>0.0 ? total_bfs_time/total_cluster_time : 0.0)<<endl;
            cout<<"[PRINT TIME] average clustering time: "<<total_cluster_time/double(query_count)<<endl;
            #ifdef EVALUATION_CLUSTERING_QUALITY
                cout<<"[PRINT ARI] average ARI score: "<<total_ari/double(query_count)<<endl;
            #endif
        }
        void print_index() override
        {
            cout<<"[PRINT HASH FOR EACH VERTEX] start !!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<i<<": "<<rnk[i]<<endl;
            }
            cout<<"[PRINT HASH FOR EACH VERTEX] finish !!!"<<endl;
            cout<<"[PRINT SKETCH] start !!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<i<<": ";
                for(auto &s:sketch[i])
                {
                    cout<<s<<' ';
                }
                cout<<endl;
            }
            cout<<"[PRINT SKETCH] finish !!!"<<endl;
            cout<<"[PRINT NEIGHBOR ORDER] start !!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<i<<": ";
                for(auto &item:neighbor_order[i])
                {
                    cout<<"("<<(float)item.first/float(PRECISION)<<","<<item.second<<")";
                }
                cout<<endl;
            }
            cout<<"[PRINT NEIGHBOR ORDER] finish !!!"<<endl;
            cout<<"[PRINT CLUSTER INDEX] start !!!"<<endl;
            for(int i=0;i<DELTA;i++)
            {
                float num=DELTA;
                cout<<"simialrity bucket ["<<i/num<<","<<(i+1)/num<<"): ";
                for(auto &item:bucket[i])
                {
                    cout<<"("<<item.first<<","<<item.second<<")";
                }
                cout<<endl;
            }
            cout<<"[PRINT CLUSTER INDEX] finish !!!"<<endl;
        }
        void store_index() override
        {
            cout<<"[STORE INDEX] start !!!!"<<endl;
            string store_no="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            string store_ci="./Index/"+dataset+"/BOTBIN_Cluster_Index.txt";
            ofstream write_index;
            write_index.open(store_no);
            for(auto &v:neighbor_order)
            {
                write_index<<v.size()<<' ';
                for(auto &item:v)
                {
                    write_index<<item.first<<' '<<item.second<<' ';
                }
                write_index<<endl;
            }
            write_index.close();
            write_index.open(store_ci);
            for(auto &b:bucket)
            {
                write_index<<b.size()<<' ';
                for(auto &item:b)
                {
                    write_index<<item.first<<' '<<item.second<<' ';
                }
                write_index<<endl;
            }
            write_index.close();
            cout<<"[STORE INDEX] finish !!!!"<<endl;
        }
        void load_index() override
        {
            cout<<"[LOAD INDEX] start !!!!"<<endl;
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
            cout<<"[LOAD INDEX] finish !!!!"<<endl;
            obtain_space_ratio();
        }

        void update(string update_type, int u, int v, int update_way) override
        {
            cout<<"[Error] current version of BOTBIN can not support update"<<endl;
        }
        void print_update_time() override
        {
            cout<<"[Error] current version of BOTBIN can not support update"<<endl;
        }    
};
