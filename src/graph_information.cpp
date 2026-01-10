/*
    this file is used to print the information of graphs, including:
        1. vertex number
        2. edge number
        3. average degree
        4. max degree
        5. average high threshold degree
        6. max high threshold degree
        7. average low threshold degree
        8. max low threshold degree
        9. the average max bucket between neighbors
        10. \overline{|PPT(v)|}
        11. |PPT|
*/
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
#include <index/PPT.h>
bool cmp_adj(int x, int y)
{
	return x>y;
}
bool cmp_order_space_index(pair<float,int> x,pair<float,int> y)
{
	return x.first>y.first;
}
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
	{"webbase","/home/hnu/Disk0/ParDon/graph/webbase.txt"},
	// {"uk2005","/home/hnu/Disk0/ParDon/graph/uk2005.txt"},
	// real graphs with billon-scale undirected edges
	// {"it2004","/home/hnu/Disk0/ParDon/graph/it2004.txt"},
	// {"friendster","/home/hnu/Disk0/ParDon/graph/friendster.txt"},
	// {"sk2005","/home/hnu/Disk0/ParDon/graph/sk2005.txt"}
};
#define INF -1
class GraphInformation
{
    private:
        string filepath;
        string dataset;
        int k=0;
        float RHO=0.1;
        float FAILURE_PB=0.001;
        int DELTA=10;
        string format;

        int** graph;
        int* degree;
        int** ht_degree;
        int** lt_degree;
        int n = 0;
        int m = 0;
        int dmax = 0;

        static vector<int> rnk;
        vector<vector<int> > sketch;

        vector<vector<pair<float, int>>> neighbor_order;
        vector<vector<int>> interval_num;
        int** max_mu_for_interval;

        vector<vector<int>> vertex2ppt;
        unordered_map<pair<int,int>,int,PairHash> point2ppt;

        vector<PPT> ppts;
        int ppt_max_id=0;
    public:
        GraphInformation(string file_path, string dataset_name)
        {
            this->filepath=file_path;
            this->dataset=dataset_name;
            read_graph();
            for(int i=0;i<n;i++)
            {
                sort(graph[i],graph[i]+degree[i],cmp_adj);
            }
            k=(int)(ceil((1/(2*pow(RHO,2)))*(log(2/FAILURE_PB))));
            construct_bottom_k_sketch();
            construct_neighbor_order();
            construct_ppt_vertex_map();
        }
        void read_graph()
        {
            cout << "[LOAD GRAPH] start!!! graph: "<<dataset<<" filepath: "<<filepath << endl;
            ifstream infile(filepath);
            if (!infile) {
                cout << "无法打开输入文件" << endl;
                return;
            }
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
            infile.close();
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
        float compute_approximate_jaccard_similarity(int u, int v)
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
            return j_sim;
        }
        void construct_bottom_k_sketch()
        {
            cout<<"[CONSTRUCT INDEX] start building BOTTOM-K-SKETCH"<<endl;
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
            cout<<"[CONSTRUCT INDEX] finish building BOTTOM-K-SKETCH"<<endl;
        }
        void construct_neighbor_order()  
        {
            cout<<"[CONSTRUCT INDEX] start building Neighbor Order"<<endl;
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
                        float similarity=compute_approximate_jaccard_similarity(u,v);
                        neighbor_order[u].push_back(make_pair(similarity,v));
                        neighbor_order[v].push_back(make_pair(similarity,u));
                        int bucket=floor(similarity*DELTA);
                        if(bucket>=DELTA) bucket=DELTA-1;
                        if(bucket < 0) bucket = 0;
                        interval_num[u][bucket]+=1;
                        interval_num[v][bucket]+=1;
                    }
                }
                sort(neighbor_order[u].begin(),neighbor_order[u].end(),cmp_order_space_index);
            }
            cout<<"[CONSTRUCT INDEX] finish building Neighbor Order"<<endl;
        }
        void construct_ppt_vertex_map()
        {
            cout<<"[CONSTRUCT Index] start PPT decomposition"<<endl;
            vertex2ppt.resize(n);
            max_mu_for_interval=new int*[n];
            for(int i=0;i<n;i++)
            {
                max_mu_for_interval[i]=new int[DELTA];
                int count=1;
                for(int bucket=DELTA-1;bucket>0;--bucket)
                {
                    count+=interval_num[i][bucket];
                    max_mu_for_interval[i][bucket]=count;
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
            cout<<"[CONSTRUCT INDEX] finish PPT decomposition"<<endl;
        }
        void print_graph_information()
        {
            cout<<"[GRAPH INFORMATION] for dataset: "<<dataset<<endl;
            cout<<"1. Vertex Number: "<<n<<endl;
            cout<<"2. Edge Number: "<<m<<endl;
            double avg_degree=double(2*m)/n;
            cout<<"3. Average Degree: "<<avg_degree<<endl;
            int max_degree=0;
            for(int i=0;i<n;i++)
            {
                if(degree[i]>max_degree)
                {
                    max_degree=degree[i];
                }
            }
            cout<<"4. Max Degree: "<<max_degree<<endl;
            ht_degree=new int*[DELTA];
            lt_degree=new int*[DELTA];
            for(int i=0;i<DELTA;i++)
            {
                ht_degree[i]=new int[n];
                lt_degree[i]=new int[n];
                memset(ht_degree[i],0,sizeof(int)*n);
                memset(lt_degree[i],0,sizeof(int)*n);
            }
            for(int i=0;i<n;i++)
            {
                for(auto item:neighbor_order[i])
                {
                    int j=item.second;
                    int neighbor_bucket=item.first*DELTA;
                    if (neighbor_bucket < 0) neighbor_bucket = 0;
                    if(neighbor_bucket>=DELTA) neighbor_bucket=DELTA-1;
                    for(int bucket=neighbor_bucket;bucket>=0;--bucket)
                    {
                        if(max_mu_for_interval[i][bucket]<max_mu_for_interval[j][bucket])
                        {
                            ht_degree[bucket][i]++;
                        }
                        if(max_mu_for_interval[i][bucket]>max_mu_for_interval[j][bucket])
                        {
                            lt_degree[bucket][i]++;
                        }
                    }
                }
            }
            int total_max_htd=0;
            int total_max_ltd=0;
            double total_avg_htd=0.0;
            double total_avg_ltd=0.0;
            for(int bucket=0;bucket<DELTA;bucket++)
            {
                double sum_htd=0.0;
                int max_htd=0;
                double sum_ltd=0.0;
                int max_ltd=0;
                for(int i=0;i<n;i++)
                {
                    sum_htd+=ht_degree[bucket][i];
                    if(ht_degree[bucket][i]>max_htd)
                    {
                        max_htd=ht_degree[bucket][i];
                    }
                    sum_ltd+=lt_degree[bucket][i];
                    if(lt_degree[bucket][i]>max_ltd)
                    {
                        max_ltd=lt_degree[bucket][i];
                    }
                }
                if(max_htd>total_max_htd)
                {
                    total_max_htd=max_htd;
                }
                if(max_ltd>total_max_ltd)
                {
                    total_max_ltd=max_ltd;
                }
                if(sum_htd != sum_ltd) cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!warning: sum_htd != sum_ltd!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                total_avg_htd+=sum_htd/n;
                total_avg_ltd+=sum_ltd/n;
                cout<<"5. Average High Threshold Degree for interval ["<<float(bucket/float(DELTA))<<","<<float((bucket+1)/float(DELTA))<<"): "<<double(sum_htd)/n<<endl;
                cout<<"6. Max High Threshold Degree for interval ["<<float(bucket/float(DELTA))<<","<<float((bucket+1)/float(DELTA))<<"): "<<max_htd<<endl;
                cout<<"7. Average Low Threshold Degree for interval ["<<float(bucket/float(DELTA))<<","<<float((bucket+1)/float(DELTA))<<"): "<<double(sum_ltd)/n<<endl;
                cout<<"8. Max Low Threshold Degree for interval ["<<float(bucket/float(DELTA))<<","<<float((bucket+1)/float(DELTA))<<"): "<<max_ltd<<endl;
            }
            total_avg_htd=total_avg_htd/DELTA;
            total_avg_ltd=total_avg_ltd/DELTA;
            cout<<"5. Overall Average High Threshold Degree: "<<total_avg_htd<<endl;
            cout<<"6. Overall Max High Threshold Degree: "<<total_max_htd<<endl;
            cout<<"7. Overall Average Low Threshold Degree: "<<total_avg_ltd<<endl;
            cout<<"8. Overall Max Low Threshold Degree: "<<total_max_ltd<<endl;
            //9. the average max bucket between neighbors
            double total_avg_max_bucket=0.0;
            for(int i=0;i<n;i++)
            {
                int max_bucket=-1;
                for(int bucket=DELTA-1;bucket>=0;--bucket)
                {
                    if(interval_num[i][bucket]>0)
                    {
                        max_bucket=bucket;
                        break;
                    }
                }
                total_avg_max_bucket+=max_bucket;
            }
            total_avg_max_bucket=total_avg_max_bucket/n;
            cout<<"9. The average max bucket between neighbors: "<<total_avg_max_bucket<<endl;
            //10. \overline{|PPT(v)|} 
            double total_avg_ppt_size=0;
            for(int i=0;i<n;i++)
            {
                total_avg_ppt_size+=vertex2ppt[i].size();
            }
            total_avg_ppt_size=total_avg_ppt_size/n;
            cout<<"10. The average |PPT(v)|: "<<total_avg_ppt_size<<endl;
            //11. |PPT|
            cout<<"11. The number of PPT: "<<ppts.size()<<endl;
        }
};
vector<int> GraphInformation::rnk=vector<int>();
int main()
{
    for(auto item: dataset_2_SCAN_file)
    {
        string dataset_name=item.first;
        string file_path=item.second;
        cout<<"================================================================"<<endl;
        cout<<"Processing dataset: "<<dataset_name<<endl;
        GraphInformation g=GraphInformation(file_path, dataset_name);
        g.print_graph_information();
    }
}