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
using namespace std;
#define UNCOMPARABLE 0xffffffff
#define EQUAL 0
#define BIG 1
#define SMALL -1
// struct Componet
// {
//     vector<int> cores;                          //core vertex
//     vector<int> non_cores;                      //non-core vertex
//     vector<int> connected_componet_ids;         //componets which must connect
//     vector<pair<int,int>> edges;                //edges corresponding to skyline connectivity
//     vector<int> connectivitys;                   //skyline connectivity
// };
// struct PPT
// {
//     int id;
//     int epsilon;
//     int mu;
//     vector<int> vertices;
//     vector<vector<int>> componet_cores;
//     vector<vector<int>> componet_non_cores;
//     vector<pair<int,int>> edges;
//     unordered_map<int,int> PPT2edges;
//     vector<int> connectivitys;
// };
struct PPT
{
    int id;
    int epsilon;
    int mu;
    vector<int> vertices;
    vector<int> componet_core_ids;
    vector<pair<int,int>> edges;
    unordered_map<int,int> cs2edges;
    vector<int> connectivitys;
};
/*
    compare PPT a and PPT b
    return value:
        0: a == b
        -1: a < b if a.epsilon < b.epsilon and a.mu < b.mu
        1: a > b if a.epsilon > b.epsilon and a.mu > b.mu
        INF: uncomparable 
*/
inline double toMB(size_t bytes) {
    return bytes / 1024.0 / 1024.0;
}
int compare_PPT(PPT & a, PPT & b)
{
    if(a.epsilon == b.epsilon && a.mu == b.mu) return EQUAL;
    if(a.epsilon > b.epsilon && a.mu > b.mu) return BIG;
    else if(a.epsilon < b.epsilon && a.mu < b.mu) return SMALL;
    else return UNCOMPARABLE;
}
int compare_similarity_SC(int & a, int & b)
{
    if(a == b ) return EQUAL;
    if(a > b) return BIG;
    else if(a < b) return SMALL;
    else return UNCOMPARABLE;
}
bool cmp_ppt(PPT & a, PPT & b)
{
    if(a.epsilon>b.epsilon)
    {
        return true;
    }
    else if(a.epsilon==b.epsilon)
    {
        if(a.mu>b.mu) return true;
        else return false;
    }
    else
    {
        return false;
    }
}
bool cmp_pair(pair<int,int> & a, pair<int,int> & b)
{
    if(a.first>b.first)
    {
        return true;
    }
    else if(a.first==b.first)
    {
        if(a.second>b.second) return true;
        else return false;
    }
    else
    {
        return false;
    }
}
struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};
namespace std {
    template<>
    struct hash<PPT> {
        size_t operator()(const PPT& p) const noexcept {
            return std::hash<int>()(p.epsilon) ^ (std::hash<int>()(p.mu) << 1);
        }
    };
}
