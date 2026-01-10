#pragma once

#include <iostream>
#include <vector>
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

class SCAN_DSU_Node
{
    public:
        int id;                 //vertex's id
        int rank;               //merge by rank
        int father;             //locate the root of componet in current level
        int hook;               /*
                                    locate subtree using vertex with the minmize mu
                                    thus, we can ensure the hook of this componet is the top node of the subtree
                                */
        int componet;           //record the higher-level's componet connectivity
        int pre_mu;             //previos mu
        int cur_mu;             //current mu
        SCAN_DSU_Node(int id)
        {
            this->id=id;
            this->rank=0;
            this->father=id;
            this->hook=id;
            this->componet=id;
            this->cur_mu=1;
        }
        void init()
        {
            this->rank=0;
            this->father=id;
            this->hook=id;
            this->componet=id;
        }
        void save_componet_and_init()
        {
            this->rank=0;
            this->father=id;
            this->hook=id;
        }
};

class SCAN_DSU
{
    public:
        vector<SCAN_DSU_Node> nodes;
        int size;
        SCAN_DSU()
        {

        }
        SCAN_DSU(int size)
        {
            this->size=size;
            this->nodes.reserve(size);
            for(int i=0;i<size;i++)
            {
                this->nodes.push_back(SCAN_DSU_Node(i));
            }
        }
        unordered_map<int,vector<int>> mu_decomposition_and_update(vector<vector<pair<float, int>>> & neighbor_order, double epsilon)
        {
            unordered_map<int,vector<int>> Vset;
            for(int i=0;i<size;i++)
            {
                int count=1;
                for(auto nei:neighbor_order[i])
                {
                    if(nei.first>=epsilon) count++;
                    else break;
                }
                nodes[i].pre_mu=nodes[i].cur_mu;
                nodes[i].cur_mu=count;
                if(Vset.find(count)!=Vset.end())
                {
                    Vset[count].push_back(i);
                }
                else
                {
                    vector<int> temp;
                    temp.push_back(i);
                    Vset[count]=temp;
                }
            }
            return Vset;
        }
        void make_set(int v)
        {
            nodes[v].init();
        }
        int find(SCAN_DSU_Node v)
        {
            if(v.father!=v.id)
            {
                //path compression
                v.father=find(nodes[v.father]);
            }
            return v.father;
        }
        void unite(int u, int v)
        {
            int ru=find(nodes[u]);
            int rv=find(nodes[v]);
            if(ru!=rv)
            {
                //merge by rank
                int componet_ru=nodes[ru].componet;
                int componet_rv=nodes[rv].componet;
                if(nodes[ru].rank<nodes[rv].rank)
                {
                    nodes[ru].father=rv;
                    /*
                        when unite u and v, this mean its tree is connected
                        ru and rv is connected
                        so they are belongs to the same mu support
                        update the new root's componet to ensure the higher mu support
                    */ 
                    if(nodes[componet_rv].cur_mu<nodes[componet_ru].cur_mu)
                    {
                        nodes[rv].componet=nodes[ru].componet;
                    }
                }
                else if(nodes[ru].rank>nodes[rv].rank)
                {
                    nodes[rv].father=ru;
                    if(nodes[componet_ru].cur_mu<nodes[componet_rv].cur_mu)
                    {
                        nodes[ru].componet=nodes[rv].componet;
                    }
                }
                else
                {
                    nodes[rv].father=ru;
                    nodes[ru].rank+=1;
                    if(nodes[componet_ru].cur_mu<nodes[componet_rv].cur_mu)
                    {
                        nodes[ru].componet=nodes[rv].componet;
                    }
                }
                // if(nodes[componet_ru].cur_mu<nodes[componet_rv].cur_mu)
                // {
                //     nodes[ru].componet=nodes[rv].componet;
                // }
                // else
                // {
                //     nodes[rv].componet=nodes[ru].componet;
                // }
            }
        }
        void update_dsu(vector<int> & V_)
        {
            for(auto v:V_)
            {
                /*
                    check the element in V
                    sync the v's componet and its root's componet
                    update the hook of the tree
                */
                int rv=find(nodes[v]);
                nodes[v].componet=nodes[rv].componet;
                int father_tree=nodes[rv].hook;
                if(nodes[father_tree].cur_mu>nodes[v].cur_mu)
                {
                    nodes[rv].hook=v;
                }
            }
        }
};
