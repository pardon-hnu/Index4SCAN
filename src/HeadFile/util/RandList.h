#pragma once
#include <algorithm>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <set>
#include <map>
using namespace std;

class RandList
{
private:
	int cap;

public:
	int vnum;
	int *vlist;
	unsigned int *vpos;
	RandList()
	{
		vlist = nullptr;
        vpos = nullptr;
	};
	void init(int _cap)
	{
		cap = _cap;
		if (vlist == nullptr)
			vlist = new int[cap];
		if (vpos == nullptr)
			vpos = new unsigned int[cap];

		vpos = new unsigned int[cap];
		vnum = 0;
		for (int i = 0; i < cap; i++)
		{
			vpos[i] = cap;
		}
	}
	inline void add(int vid)
	{
		vlist[vnum] = vid;
		vpos[vid] = vnum++;
	};
	inline void remove(int vid)
	{
		int last_id = vlist[vnum - 1];
		unsigned int id_pos = vpos[vid];
		vlist[id_pos] = last_id;
		vpos[last_id] = id_pos;

		--vnum;
		vpos[vid] = cap;
	}
	void clear()
	{
		for (int i = 0; i < vnum; i++)
			vpos[vlist[i]] = cap;
		vnum = 0;
	}
	int get(int i)
	{
		return vlist[i];
	}
	bool contains(int vid)
	{
		return vpos[vid] != cap;
	}
	bool empty() { return vnum == 0; }
	int getSize() { return vnum; }
	int getCap() { return cap; }
	void dispose()
	{
		if (vlist != nullptr)
		{
			delete[] vlist;
			vlist = nullptr;
		}
		if (vpos != nullptr)
		{
			delete[] vpos;
			vpos = nullptr;
		}
	}
	~RandList()
	{
		dispose();
	}
};


class UnionFind {
public:
    int node_num;
    vector<int> father;
    vector<int> rank;
    UnionFind()=default;
    UnionFind(int n);
    void Init();
    void Init(vector<int> & nums);
    void Init(RandList & nums);
    int Find(int x);
    void Unite(int x, int y);
    int UniteRoot(int x, int y);
    void BatchUnite(set<int> &nodes);
};


UnionFind::UnionFind(int n) {
    node_num = n;
    father.resize(node_num);
    rank.resize(node_num);
    for (int i = 0; i < node_num; i++) {
        father[i] = i;
        rank[i] = 0;
    }
}

void UnionFind::Init() {
    int i = 0;
    for (int& element : father) {
        element = i++;
    }
    for (int& element : rank) {
        element = 0;
    }
}
void UnionFind::Init(vector<int> &nums) {
    for (auto & element : nums) {
        father[element] = element;
        rank[element] = 0;
    }
}
void UnionFind::Init(RandList & nums){
    for(int i=0;i<nums.vnum;i++)
    {
        int element=nums.vlist[i];
        father[element] = element;
        rank[element] = 0;
    }
}

int UnionFind::Find(int x) {
    if(x == father[x])
    return x;
    return father[x] = Find(father[x]); 
}

void UnionFind::Unite(int x, int y){
    x = Find(x);
    y = Find(y);
    if(x == y)
        return ;
    if(rank[x] < rank[y]){
        father[x] = y;
    }else{
        father[y] = x;
        if(rank[x] == rank[y])
            rank[x]++;
    }
}

int UnionFind::UniteRoot(int x, int y) {
    if (x == y)
        return x;
    if (rank[x] < rank[y]) {
        father[x] = y;
        return y;
    } else {
        father[y] = x;
        if(rank[x] == rank[y])
            rank[x]++;
        return x;
    }
}

//void UnionFind
void UnionFind::BatchUnite(set<int> &root_nodes) {
    if (root_nodes.size() == 0)
        return;
    vector<int> vec_root_nodes;
    for (int node : root_nodes) {
        vec_root_nodes.push_back(node);
    }
    for (int i = 0; i < vec_root_nodes.size() - 1; i++) {
        int x = vec_root_nodes[i];
        int y = vec_root_nodes[i+1];
        if(rank[x] < rank[y]){
            father[x] = y;
        }else{
            father[y] = x;
            if(rank[x] == rank[y])
                rank[x]++;
        }
    }
}
