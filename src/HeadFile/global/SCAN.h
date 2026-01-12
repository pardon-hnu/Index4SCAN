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
#include <index/GS_Index.h>
#include <index/BOTBIN.h>
#include <index/Forest_Index.h>
#include <index/PPT_Index.h>
#include <memory>
using namespace std;
typedef std::chrono::duration<double> tms;

class SCAN
{
    private:
        string dataset;
        string method;
        std::unique_ptr<SCAN_Index> index;
    public:
        SCAN(string dataset,string method,vector<string> parameters)
        {
            this->dataset=dataset;
            this->method=method;
            if(method=="GS-Index")
            {
                index = std::make_unique<GS_Index>(dataset);
            }
            else if(method=="BOTBIN")
            {
                if(parameters.size())
                {
                    float RHO=stof(parameters[0]);
                    float FAILURE_PB=stof(parameters[1]);
                    int DELTA=stoi(parameters[2]);
                    index = std::make_unique<BOTBIN>(dataset,RHO,FAILURE_PB,DELTA);
                }
                else
                {
                    index = std::make_unique<BOTBIN>(dataset);
                }
            }
            else if(method=="FOREST")
            {
                // index = std::make_unique<SCAN_Forest_Index>(dataset);
                if(parameters.size())
                {
                    int DELTA=stoi(parameters[0]);
                    index = std::make_unique<SCAN_Forest_Index>(dataset,DELTA);
                }
                else
                {
                    index = std::make_unique<SCAN_Forest_Index>(dataset);
                }
            }
            else if(method=="PPT")
            {
                cout << "parameters.size()=" << parameters.size() << endl;
                for(auto &p: parameters) cout << "[" << p << "] ";
                cout << endl;
                if(parameters.size()==0)
                {
                    index = std::make_unique<SCAN_PPT_Index>(dataset);
                }
                else if(parameters.size()==1)
                {
                    int p=stoi(parameters[0]);
                    if(p>=-3&&p<=-1)
                    {
                        int format_id=p;
                        index = std::make_unique<SCAN_PPT_Index>(dataset,10,format_id);
                    }
                    else if(p>0)
                    {
                        int DELTA=p;
                        index = std::make_unique<SCAN_PPT_Index>(dataset,DELTA,-2);
                    }
                    else
                    {
                        cout<<"[Error] parmeters size: 1 and parameter error !!!"<<endl;
                    }
                }
                else if(parameters.size()==2)
                {
                    int DELTA=stoi(parameters[0]);
                    int format_id=stoi(parameters[1]);
                    if(format_id<-3||format_id>-1)
                    {
                        cout<<"[Error] parmeters size: 2 and parmeters error !!!"<<endl;
                    }
                    else
                    {
                        index = std::make_unique<SCAN_PPT_Index>(dataset,DELTA,format_id);
                    }
                }
                else
                {
                    cout<<"[Error] parmeters size >2 and error !!!"<<endl;
                }
            }
            else
            {
                cout<<"[Error] there no such index !!!!"<<endl;
            }
        }
        ~SCAN()
        {

        }
        void index_test()
        {
            index->load_index();
            index->print_index();
            index->query(0.4,2);
        }
        void index_construct()
        {
            index->construct();
            index->print_build_cost();
            index->store_index();
        }
        void single_cluster(float epsilon,int mu)
        {
            index->load_index();
            index->query(epsilon,mu);
            string output_cluster_path="Result/result/"+dataset+"/"+method+"-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".txt";
            index->write_result(output_cluster_path);
        }
        void batch_cluster(vector<float> epsilons,vector<int> mus, bool write_result=false, bool write_labels=false)
        {
            index->load_index();
            for(auto epsilon:epsilons)
            {
                for(auto mu:mus)
                {
                    cout<<"epsilon: "<<epsilon<<" mu: "<<mu<<endl;
                    index->query(epsilon,mu);
                    if(write_result)
                    {
                        string output_cluster_path="Result/result/"+dataset+"/"+method+"-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".txt";
                        index->write_result(output_cluster_path);
                    }
                    if(write_labels)
                    {
                        string output_labels_path="Result/labels/"+dataset+"/"+method+"-"+std::to_string(epsilon)+"-"+std::to_string(mu)+".bin";
                        index->write_binary_labels(output_labels_path);
                    }
                }
            }
            index->print_cluster_time();
        }
};