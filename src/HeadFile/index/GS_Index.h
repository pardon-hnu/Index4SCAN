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
using namespace std;
typedef std::chrono::duration<double> tms;
bool cmp_order(pair<int,int> x,pair<int,int> y)
{
	return x.second>y.second;
}
class GS_Index : public SCAN_Index{
    private:
        vector<vector<pair<int,int>>> neighbor_order;
        vector<vector<pair<int,int>>> core_order;
        int query_count=0;
    public:
        GS_Index(string dataset)
        {
            init_graph(dataset);
        }
        ~GS_Index()
        {

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
            cout<<"[CONSTRUCT INDEX] begin constructing neighbor order!!!"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            neighbor_order.resize(n);
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<degree[i];j++)
                {
                    int v=graph[i][j];
                    if(v>i)
                    {
                        int score=compute_similarity_exact(i,v);
                        neighbor_order[i].push_back(make_pair(v,score));
                        neighbor_order[v].push_back(make_pair(i,score));
                    }
                }
                sort(neighbor_order[i].begin(),neighbor_order[i].end(),cmp_order);
            }
            const auto end=std::chrono::steady_clock().now();
            tms no_construct_time=end-begin;
            cout<<"[COSNTRUCT INDEX] finishi construting no !!!"<<endl;
            cout<<"[TIME COST] "<<no_construct_time.count()<<" s."<<endl;
            total_build_time+=no_construct_time.count();
            for(auto &item:neighbor_order)
            {
                total_build_space+=item.size()*8;
            }
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
                    item.push_back(make_pair(v,sim));
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
        void construct_core_order()
        {
            cout<<"[CONSTRUCT INDEX] begin constructing core order !!!"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            core_order.resize(dmax);
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<neighbor_order[i].size();j++)
                {
                    core_order[j].push_back(make_pair(i,neighbor_order[i][j].second));
                }
            }
            for(int i=0;i<dmax;i++)
            {
                sort(core_order[i].begin(),core_order[i].end(),cmp_order);
            }
            const auto end=std::chrono::steady_clock().now();
            tms co_construct_time=end-begin;
            cout<<"[CONSTRUCT INDEX] finish constructing co !!!"<<endl;
            cout<<"[TIME COST] "<<co_construct_time.count()<<" s."<<endl;
            total_build_time+=co_construct_time.count();
            for(auto &item:core_order)
            {
                total_build_space+=item.size()*8;
            }
        }
        void clustering_gs_index(int epsilon, int mu)
        {
            clear_cluster_result();
            cout<<"[CLUSTER] clustering with GS-index start !!!"<<endl;
            const auto begin=std::chrono::steady_clock().now();
            for(auto &core:core_order[mu-2])
            {
                if(labels[core.first]!=-1) continue;
                if(core.second<epsilon) break;
                // cout<<"iterate vertex "<<core.first<<endl;
                queue<int> cluster;
                unordered_set<int> cluster_set;
                cluster.push(core.first);
                cluster_set.insert(core.first);
                while(!cluster.empty())
                {
                    int v=cluster.front();
                    cluster.pop();
                    labels[v]=cluster_number;
                    for(auto &nv:neighbor_order[v])
                    {
                        if(labels[nv.first]!=-1) continue;
                        if(nv.second<epsilon) break;
                        // cout<<"add its neighbors "<<nv.first<<" into clsuter"<<endl;
                        if(neighbor_order[nv.first].size()>=mu-1&&neighbor_order[nv.first][mu-2].second>=epsilon) //core
                        {
                            // cout<<"core "<<endl;
                            if(cluster_set.find(nv.first)==cluster_set.end())
                            {
                                cluster.push(nv.first);
                                cluster_set.insert(nv.first);
                            }
                        }
                        else //non core
                        {
                            // cout<<"non_core "<<endl;
                            labels[nv.first]=cluster_number;
                        }
                    }
                }
                cluster_number++;
            }
            const auto end=std::chrono::steady_clock().now();
            tms cluster_time=end-begin;
            cout<<"[CLSUTER] finish clustering!!!"<<endl;
            cout<<"[TIME COST] "<<cluster_time.count()<<" s."<<endl;
            total_cluster_time+=cluster_time.count();
            query_count++;
        }
        void construct() override
        {
            #ifdef FAST_CONSTRUCT
                fast_construct_neighbor_order();
            #else
                construct_neighbor_order();
            #endif
		    construct_core_order();
        }
        void query(float epsilon, int mu) override
        {
            int eps=epsilon*PRECISION;
            clustering_gs_index(eps, mu);
        }
        void print_build_cost() override
        {
            cout<<"[PRINT TIME] total building time: "<<total_build_time<<endl;
            cout<<"[PRINT SPAEC]: total space cost: "<<(total_build_space/1024)/1024<<" MB"<<endl;
        }
        void print_cluster_time() override
        {
            cout<<"[PRINT TIME] total clustering time: "<<total_cluster_time<<endl;
            cout<<"[PRINT TIME] average clustering time: "<<total_cluster_time/double(query_count)<<endl;
        }
        void print_index() override
        {
            cout<<"[PRINT NEIGHBOR ORDER] start!!!"<<endl;
            for(int i=0;i<n;i++)
            {
                cout<<"vertex "<<i<<" and its nerighbor order: ";
                for(int j=0;j<degree[i];j++)
                {
                    cout<<"("<<neighbor_order[i][j].first<<","<<float(neighbor_order[i][j].second)/float(PRECISION)<<") ";
                }
                cout<<endl;
            }
            cout<<"[PRINT NEIGHBOR ORDER] finish!!!"<<endl;
            cout<<"[PRINT CORE ORDER] start !!!"<<endl;
            for(int i=0;i<dmax;i++)
            {
                cout<<"CO["<<i+2<<"]:";
                for(int j=0;j<core_order[i].size();j++)
                {
                    cout<<"("<<core_order[i][j].first<<","<<float(core_order[i][j].second)/float(PRECISION)<<") ";
                }
                cout<<endl;
            }
            cout<<"[PRINT CORE ORDER] finish !!!"<<endl;
        }
        void store_index() override
        {
            cout<<"[STORE INDEX] start !!!!"<<endl;
            string store_no="./Index/"+dataset+"/GS_Index_Neighbor_Order.txt";
            string store_co="./Index/"+dataset+"/GS_Index_Core_Order.txt";
            ofstream write_index;
            write_index.open(store_no);
            for(auto &item:neighbor_order)
            {
                write_index<<item.size()<<' ';
                for(auto &p:item)
                {
                    write_index<<p.first<<' '<<p.second<<' ';
                }
                write_index<<endl;
            }
            write_index.close();
            write_index.open(store_co);
            for(auto &item:core_order)
            {
                write_index<<item.size()<<' ';
                for(auto &p:item)
                {
                    write_index<<p.first<<' '<<p.second<<' ';
                }
                write_index<<endl;
            }
            write_index.close();
            cout<<"[STORE INDEX] finish !!!!"<<endl;
        }
        void load_index() override
        {
            cout<<"[LOAD INDEX] start loading gs index !!!!"<<endl;
            string load_no="./Index/"+dataset+"/GS_Index_Neighbor_Order.txt";
            string load_co="./Index/"+dataset+"/GS_Index_Core_Order.txt";
            string load_botbin_no="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            ifstream read_index;
            string line;
            if(dataset=="sk2005") 
            {
                read_index.open(load_botbin_no);
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
                        item.push_back(make_pair(v,sim));
                    }
                    neighbor_order.push_back(item);
                }
                read_index.close();
            }
            else
            {
                read_index.open(load_no);
                while(getline(read_index, line))
                {
                    char* temp;
                    vector<pair<int,int>> item;
                    temp = strtok(const_cast<char*>(line.c_str()), " ");
                    int len=atoi(temp);
                    for(int i=0;i<len;i++)
                    {
                        temp = strtok(NULL, " ");
                        int v=atoi(temp);
                        temp = strtok(NULL, " ");
                        int sim=atoi(temp);
                        item.push_back(make_pair(v,sim));
                    }
                    neighbor_order.push_back(item);
                }
                read_index.close();
            }
            read_index.open(load_co);
            while(getline(read_index, line))
            {
                char* temp;
                vector<pair<int,int>> item;
                temp = strtok(const_cast<char*>(line.c_str()), " ");
                int len=atoi(temp);
                for(int i=0;i<len;i++)
                {
                    temp = strtok(NULL, " ");
                    int v=atoi(temp);
                    temp = strtok(NULL, " ");
                    int sim=atoi(temp);
                    item.push_back(make_pair(v,sim));
                }
                core_order.push_back(item);
            }
            read_index.close();
            cout<<"[LOAD INDEX] finish loading gs index !!!!"<<endl;
        }
};