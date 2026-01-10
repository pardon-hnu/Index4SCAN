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
typedef std::chrono::duration<double> tms;
bool cmp_adj(int x, int y)
{
	return x>y;
}
class SCAN_Index
{
    public:
        //graph information
        int** graph;
        int* degree;
        int n = 0;
        int m = 0;
        int dmax = 0;
        string dataset;
        //clusters result
        vector<int> labels;
        int cluster_number=0;
        double total_cluster_time=0.0;
        double total_build_time=0.0;
        double total_build_space=0.0;
        //shared function
        void init_graph(string dataset)
        {   
            this->dataset = dataset;
            string filepath = dataset_2_SCAN_file[dataset];
            read_graph(filepath);
            sort_adj();
            for(int i=0;i<n;i++)
            {
                if(degree[i]>dmax) dmax=degree[i];
            }
            labels.reserve(n);
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
        void clear_cluster_result()
        {
            labels.clear();
            for(int i=0;i<n;i++)
            {
                labels.push_back(INF);
            }
            cluster_number=0;
        }
        void print_graph()
        {
            cout << "[PRINT GRAPH] start!!!" << endl;
            cout << "graph has " << n << " vertices and " << m << " edges" << endl;
            for (int i = 0;i < n;i++)
            {
                cout << "id: " << i << " deg: " << degree[i] << " neighbors: ";
                for (int j = 0;j < degree[i];j++)
                {
                    cout << graph[i][j] << ' ';
                }
                cout << endl;
            }
            cout << "[PRINT GRAPH] success!!!" << endl;
        }
        void print_result()
        {
            cout<<"[PRINT CLSUTER RESUTL] start !!!"<<endl;
            map<int,vector<int>> clusters;
            for(int i=0;i<n;i++)
            {
                clusters[labels[i]].push_back(i);
            }
            for(auto &cluster:clusters)
            {
                cout<<"cluster ("<<cluster.first<<"): ";
                for(auto &v:cluster.second)
                {
                    cout<<v<<' ';
                }
                cout<<endl;
            }
            cout<<"[PRINT CLSUTER RESUTL] finish !!!"<<endl;
        }
        void write_labels(string filename)
        {
            cout<<"[WRITE LABELS] start !!!"<<endl;
            ofstream write_label;
            write_label.open(filename);
            for(int i=0;i<n;i++)
            {
                write_label<<labels[i]<<' ';
            }
            write_label.close();
            cout<<"[WRITE LABELS] finish !!!"<<endl;
        }
        void write_binary_labels(string filename)
        {
            cout<<"[WRITE LABELS] start !!!"<<endl;
            ofstream write_label(filename, ios::binary);
            if (!write_label) 
            {
                cerr << "Cannot open file: " << filename << endl;
                exit(1);
            }
            write_label.write(reinterpret_cast<const char*>(labels.data()),sizeof(int)*n);
            write_label.close();
            cout<<"[WRITE LABELS] finish !!!"<<endl;
        }
        void write_result(string filename)
        {
            cout<<"[WRITE CLSUTER RESUTL] start !!!"<<endl;
            map<int,vector<int>> clusters;
            for(int i=0;i<n;i++)
            {
                clusters[labels[i]].push_back(i);
            }
            ofstream write_cluster;
            write_cluster.open(filename);
            for(auto &cluster:clusters)
            {
                if(cluster.first==-1)
                {
                    continue;
                }
                for(auto &v:cluster.second)
                {
                    write_cluster<<v<<' ';
                }
                write_cluster<<endl;
            }
            write_cluster.close();
            cout<<"[WRITE CLSUTER RESUTL] finish !!!"<<endl;
        }
        virtual void construct()=0;
        virtual void query(float epsilon, int mu)=0;
        virtual void print_index()=0;
        virtual void print_build_cost()=0;
        virtual void print_cluster_time()=0;
        virtual void load_index()=0;
        virtual void store_index()=0;
        virtual ~SCAN_Index() = default;
};

class DSCAN_Index
{

};