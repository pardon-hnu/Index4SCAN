/*
    Dataset statistics for the space and construction complexity analyses.
    Reuse stored similarities; do not reconstruct sketches or indexes.
    Usage: ./complexity_information [dataset|all] [delta]
*/
#include "HeadFile/global/global.h"

class ComplexityInformation
{
    private:
        string dataset;
        int delta;
        int n=0;
        long long m=0;
        vector<int> ct;
        vector<long long> tree_vertices;
        vector<long long> tree_edges;
        vector<long long> unequal_edges;

    public:
        ComplexityInformation(string dataset_name, int interval_num)
        {
            dataset=dataset_name;
            delta=interval_num;
        }

        bool collect_graph_information()
        {
            string filepath=dataset_2_SCAN_file.at(dataset);
            ifstream read_graph(filepath);
            if(!(read_graph>>n>>m))
            {
                cerr<<"[LOAD ERROR] "<<filepath<<endl;
                return false;
            }
            read_graph.close();

            string order_path="./Index/"+dataset+"/BOTBIN_Neighbor_Order.txt";
            ifstream read_order(order_path);
            if(!read_order)
            {
                cerr<<"[LOAD ERROR] "<<order_path<<endl;
                return false;
            }
            cout<<"[LOAD GRAPH] dataset: "<<dataset<<" filepath: "<<filepath<<endl;
            cout<<"[LOAD INDEX] "<<order_path<<endl;
            ct.assign(size_t(n)*delta,0);
            tree_vertices.assign(delta,0);
            tree_edges.assign(delta,0);
            unequal_edges.assign(delta,0);
            int bucket_width=PRECISION/delta;

            // First pass: histogram each row, then recover closed-neighborhood CT.
            // Include bucket 0 even though the serialized Forest only has a placeholder.
            for(int u=0;u<n;u++)
            {
                int degree;
                read_order>>degree;
                size_t offset=size_t(u)*delta;
                for(int j=0;j<degree;j++)
                {
                    int similarity,v;
                    read_order>>similarity>>v;
                    int bucket=min(similarity/bucket_width,delta-1);
                    ct[offset+bucket]++;
                }
                int count=1;
                for(int bucket=delta-1;bucket>=0;bucket--)
                {
                    count+=ct[offset+bucket];
                    ct[offset+bucket]=count;
                    if(count>=2)
                    {
                        tree_vertices[bucket]++;
                    }
                }
            }
            if(!read_order)
            {
                cerr<<"[LOAD ERROR] incomplete Neighbor Order: "<<order_path<<endl;
                return false;
            }
            read_order.close();

            // Second pass: count each undirected edge once, using final CT at both ends.
            // No adjacency lists or O(m) edge buffer are retained in memory.
            read_order.open(order_path);
            long long stored_edges=0;
            for(int u=0;u<n;u++)
            {
                int degree;
                read_order>>degree;
                for(int j=0;j<degree;j++)
                {
                    int similarity,v;
                    read_order>>similarity>>v;
                    if(u>=v) continue;
                    stored_edges++;
                    int max_bucket=min(similarity/bucket_width,delta-1);
                    for(int bucket=0;bucket<=max_bucket;bucket++)
                    {
                        tree_edges[bucket]++;
                        if(ct[size_t(u)*delta+bucket]!=ct[size_t(v)*delta+bucket])
                        {
                            unequal_edges[bucket]++;
                        }
                    }
                }
            }
            if(!read_order)
            {
                cerr<<"[LOAD ERROR] incomplete Neighbor Order: "<<order_path<<endl;
                return false;
            }
            if(stored_edges!=m)
            {
                cerr<<"[WARNING] graph header edges: "<<m
                    <<" Neighbor Order edges: "<<stored_edges<<endl;
            }
            vector<int>().swap(ct);
            return true;
        }

        void print_graph_information()
        {
            cout<<fixed<<setprecision(6);
            cout<<"[GRAPH INFORMATION] for dataset: "<<dataset<<endl;
            cout<<"1. Vertex Number: "<<n<<endl;
            cout<<"2. Edge Number: "<<m<<endl;
            cout<<"3. Delta: "<<delta<<endl;
            cout<<"[PER TREE] tree_id epsilon n_tree m_tree m_unequal_tree"<<endl;
            for(int bucket=0;bucket<delta;bucket++)
            {
                cout<<bucket<<" "<<double(bucket)/delta<<" "<<tree_vertices[bucket]
                    <<" "<<tree_edges[bucket]<<" "<<unequal_edges[bucket]<<endl;
            }

            // Report both conventions explicitly; never average the 0 0 placeholder.
            for(int first_bucket=0;first_bucket<=1;first_bucket++)
            {
                long long vertex_sum=0,edge_sum=0,unequal_sum=0;
                for(int bucket=first_bucket;bucket<delta;bucket++)
                {
                    vertex_sum+=tree_vertices[bucket];
                    edge_sum+=tree_edges[bucket];
                    unequal_sum+=unequal_edges[bucket];
                }
                int tree_num=delta-first_bucket;
                cout<<(first_bucket==0 ? "[PAPER: ALL BUCKETS]" : "[INDEX: EXCLUDE TREE 0]")<<endl;
                cout<<"4. n_tree_avg: "<<double(vertex_sum)/tree_num<<endl;
                cout<<"5. m_tree_avg: "<<double(edge_sum)/tree_num<<endl;
                cout<<"6. m_unequal_tree_avg: "<<double(unequal_sum)/tree_num<<endl;
                cout<<"7. n_interval_avg(V): "<<(n ? double(vertex_sum)/n : 0)<<endl;
                cout<<"8. Sum of vertex intervals: "<<vertex_sum<<endl;
            }
        }

        void print_index_information()
        {
            // Stream serialized records without materializing PPTs, slices, or trees.
            string ppt_path="./Index/"+dataset+"/ppts_"+to_string(delta)+".txt";
            ifstream read_ppt(ppt_path);
            if(read_ppt)
            {
                long long ppt_num,cs_num,con_num=0,core_num=0,nc_num=0;
                read_ppt>>ppt_num>>cs_num;
                for(long long i=0;i<ppt_num;i++)
                {
                    int id,epsilon,mu;
                    long long slice_num,edge_num;
                    read_ppt>>id>>epsilon>>mu>>slice_num>>edge_num;
                    con_num+=edge_num;
                    int value;
                    for(long long j=0;j<slice_num+3*edge_num;j++)
                    {
                        read_ppt>>value;
                    }
                }
                for(int list_type=0;list_type<2;list_type++)
                {
                    for(long long i=0;i<cs_num;i++)
                    {
                        long long len;
                        read_ppt>>len;
                        if(list_type==0) core_num+=len;
                        else nc_num+=len;
                        int value;
                        for(long long j=0;j<len;j++) read_ppt>>value;
                    }
                }
                if(!read_ppt)
                {
                    cerr<<"[LOAD ERROR] incomplete PPT index: "<<ppt_path<<endl;
                }
                else
                {
                    long long interval_sum=accumulate(tree_vertices.begin(),tree_vertices.end(),0LL);
                    cout<<"[PPT INFORMATION] "<<ppt_path<<endl;
                    cout<<"9. PPT Entries: "<<ppt_num<<endl;
                    cout<<"10. Cluster Slices: "<<cs_num<<endl;
                    cout<<"11. Core Occurrences (sum |PPT(v)| in stored index): "<<core_num<<endl;
                    cout<<"12. Non-core Occurrences: "<<nc_num<<endl;
                    cout<<"13. Con Records: "<<con_num<<endl;
                    cout<<"14. Con / Paper Interval Bound: "
                        <<(interval_sum ? double(con_num)/interval_sum : 0)<<endl;
                }
            }
            else cout<<"[SKIP] PPT index not found: "<<ppt_path<<endl;

            string forest_path="./Index/"+dataset+"/FOREST_Forest_"+to_string(delta)+".txt";
            ifstream read_forest(forest_path);
            if(!read_forest)
            {
                cout<<"[SKIP] Forest index not found: "<<forest_path<<endl;
                return;
            }
            cout<<"[FOREST INFORMATION] tree_id nodes_including_dummy core_records nc_records"<<endl;
            int tree_id;
            long long node_num;
            while(read_forest>>tree_id>>node_num)
            {
                long long core_num=0,nc_num=0;
                for(long long i=0;i<node_num;i++)
                {
                    int id,mu,father,value;
                    read_forest>>id>>mu>>father;
                    for(int list_type=0;list_type<3;list_type++)
                    {
                        long long len;
                        read_forest>>len;
                        if(list_type==1) core_num+=len;
                        if(list_type==2) nc_num+=len;
                        for(long long j=0;j<len;j++) read_forest>>value;
                    }
                }
                if(!read_forest)
                {
                    cerr<<"[LOAD ERROR] incomplete Forest index: "<<forest_path<<endl;
                    return;
                }
                cout<<tree_id<<" "<<node_num<<" "<<core_num<<" "<<nc_num<<endl;
                if(node_num && core_num!=tree_vertices[tree_id])
                {
                    cout<<"[WARNING] Forest core count differs from Neighbor Order at tree "<<tree_id<<endl;
                }
            }
        }
};

int main(int argc, char** argv)
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string selected_dataset=argc>1 ? argv[1] : "all";
    int delta=argc>2 ? atoi(argv[2]) : 10;
    if(delta<2 || delta>PRECISION || PRECISION%delta!=0)
    {
        cerr<<"Delta must divide PRECISION and be between 2 and "<<PRECISION<<endl;
        return 1;
    }
    if(selected_dataset!="all" && !dataset_2_SCAN_file.count(selected_dataset))
    {
        cerr<<"Unknown dataset: "<<selected_dataset<<endl;
        return 1;
    }
    int failed=0;
    for(auto &item:dataset_2_SCAN_file)
    {
        if(selected_dataset!="all" && selected_dataset!=item.first) continue;
        cout<<"================================================================"<<endl;
        cout<<"Processing dataset: "<<item.first<<endl;
        ComplexityInformation information(item.first,delta);
        if(!information.collect_graph_information())
        {
            failed++;
            continue;
        }
        information.print_graph_information();
        information.print_index_information();
    }
    return failed ? 1 : 0;
}
