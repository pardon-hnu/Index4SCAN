#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <map>
using namespace std;
#define PRECISION 100
#define INF -1
#define FAST_CONSTRUCT
#define EVALUATION_CLUSTERING_QUALITY
map<string, string> dataset_2_SCAN_file = {
	//tiny graphs
	{"example","../datasets/example.txt"},
	{"hiv","../datasets/hiv.txt"},
	{"dolphins","../datasets/dolphins.txt"},
	{"contiguous","../datasets/contiguous.txt"},
	{"club","../datasets/zachary_karate_club.txt"},
	{"tribes","/home/hnu/Disk0/ParDon/graph/Tribes.txt"},
	//real graphs
	{"amazon","/home/hnu/Disk0/ParDon/graph/amazon.txt"},
	{"dblp","/home/hnu/Disk0/ParDon/graph/dblp.txt"},
	{"youtube","/home/hnu/Disk0/ParDon/graph/youtube.txt"},
	{"skitter","/home/hnu/Disk0/ParDon/graph/skitter.txt"},
	{"pokec","/home/hnu/Disk0/ParDon/graph/pokec.txt"},
	{"topcats","/home/hnu/Disk0/ParDon/graph/topcats.txt"},
	{"livejournal1","/home/hnu/Disk0/ParDon/graph/livejournal1.txt"},
	{"livejournal2","/home/hnu/Disk0/ParDon/graph/livejournal2.txt"},
	{"orkut1","/home/hnu/Disk0/ParDon/graph/orkut1.txt"},
	{"orkut2","/home/hnu/Disk0/ParDon/graph/orkut2.txt"},
	{"indochina","/home/hnu/Disk0/ParDon/graph/indochina.txt"},
	{"uk2002","/home/hnu/Disk0/ParDon/graph/uk2002.txt"},
	//real graphs with billion-scale directed edges
	{"webbase","/home/hnu/Disk0/ParDon/graph/webbase.txt"},
	{"uk2005","/home/hnu/Disk0/ParDon/graph/uk2005.txt"},
	//real graphs with billon-scale undirected edges
	{"it2004","/home/hnu/Disk0/ParDon/graph/it2004.txt"},
	{"friendster","/home/hnu/Disk0/ParDon/graph/friendster.txt"},
	{"sk2005","/home/hnu/Disk0/ParDon/graph/sk2005.txt"}
};

static inline int64_t comb2(int64_t x) {
    return (x >= 2) ? x * (x - 1) / 2 : 0;
}

double adjusted_rand_index(const vector<int> &label1,
                           const vector<int> &label2)
{
    const size_t n = label1.size();
    if (n != label2.size() || n < 2) return 0.0;

    /*--------------------------------------------------
     * Step 1: compress labels (important if labels are large)
     *--------------------------------------------------*/
    vector<int> l1 = label1, l2 = label2;

    auto compress = [](vector<int> &v) {
        vector<int> tmp = v;
        sort(tmp.begin(), tmp.end());
        tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
        for (auto &x : v) {
            x = lower_bound(tmp.begin(), tmp.end(), x) - tmp.begin();
        }
        return (int)tmp.size();
    };

    int c1 = compress(l1);
    int c2 = compress(l2);

    /*--------------------------------------------------
     * Step 2: count a_i and b_j
     *--------------------------------------------------*/
    vector<int64_t> a(c1, 0), b(c2, 0);
    for (size_t i = 0; i < n; ++i) {
        ++a[l1[i]];
        ++b[l2[i]];
    }

    int64_t sum_ai = 0, sum_bj = 0;
    for (int i = 0; i < c1; ++i) sum_ai += comb2(a[i]);
    for (int j = 0; j < c2; ++j) sum_bj += comb2(b[j]);

    /*--------------------------------------------------
     * Step 3: count n_ij via sorting pairs
     *--------------------------------------------------*/
    vector<uint64_t> pairs;
    pairs.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        // pack (l1, l2) into one 64-bit integer
        uint64_t key = (uint64_t(l1[i]) << 32) | uint32_t(l2[i]);
        pairs.push_back(key);
    }

    sort(pairs.begin(), pairs.end());

    int64_t sum_nij = 0;
    for (size_t i = 0; i < n; ) {
        size_t j = i + 1;
        while (j < n && pairs[j] == pairs[i]) ++j;
        sum_nij += comb2(j - i);
        i = j;
    }

    /*--------------------------------------------------
     * Step 4: ARI formula
     *--------------------------------------------------*/
    const double total_pairs = double(n) * (n - 1) / 2.0;
    const double expected = (double)sum_ai * (double)sum_bj / total_pairs;
    const double max_index = 0.5 * (sum_ai + sum_bj);
    const double numerator = sum_nij - expected;
    const double denominator = max_index - expected;

    if (denominator == 0.0) return 0.0;
    return numerator / denominator;
}

double adjusted_rand_index_ignore_minus1(
    const vector<int> &label1,
    const vector<int> &label2)
{
    const size_t n0 = label1.size();
    if (n0 != label2.size() || n0 < 2) return 0.0;

    /*--------------------------------------------------
     * Step 0: filter out points with label -1
     *--------------------------------------------------*/
    vector<int> f1, f2;
    f1.reserve(n0);
    f2.reserve(n0);

    for (size_t i = 0; i < n0; ++i) {
        if (label1[i] != -1 && label2[i] != -1) {
            f1.push_back(label1[i]);
            f2.push_back(label2[i]);
        }
    }

    const size_t n = f1.size();
    if (n < 2) return 1.0;

    /*--------------------------------------------------
     * Step 1: compress labels
     *--------------------------------------------------*/
    vector<int> l1 = f1, l2 = f2;

    auto compress = [](vector<int> &v) {
        vector<int> tmp = v;
        sort(tmp.begin(), tmp.end());
        tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
        for (auto &x : v) {
            x = lower_bound(tmp.begin(), tmp.end(), x) - tmp.begin();
        }
        return (int)tmp.size();
    };

    int c1 = compress(l1);
    int c2 = compress(l2);

    /*--------------------------------------------------
     * Step 2: count a_i and b_j
     *--------------------------------------------------*/
    vector<int64_t> a(c1, 0), b(c2, 0);
    for (size_t i = 0; i < n; ++i) {
        ++a[l1[i]];
        ++b[l2[i]];
    }

    int64_t sum_ai = 0, sum_bj = 0;
    for (int i = 0; i < c1; ++i) sum_ai += comb2(a[i]);
    for (int j = 0; j < c2; ++j) sum_bj += comb2(b[j]);

    /*--------------------------------------------------
     * Step 3: count n_ij
     *--------------------------------------------------*/
    vector<uint64_t> pairs;
    pairs.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        uint64_t key = (uint64_t(l1[i]) << 32) | uint32_t(l2[i]);
        pairs.push_back(key);
    }

    sort(pairs.begin(), pairs.end());

    int64_t sum_nij = 0;
    for (size_t i = 0; i < n; ) {
        size_t j = i + 1;
        while (j < n && pairs[j] == pairs[i]) ++j;
        sum_nij += comb2(j - i);
        i = j;
    }

    /*--------------------------------------------------
     * Step 4: ARI
     *--------------------------------------------------*/
    const double total_pairs = double(n) * (n - 1) / 2.0;
    const double expected = (double)sum_ai * (double)sum_bj / total_pairs;
    const double max_index = 0.5 * (sum_ai + sum_bj);

    const double denom = max_index - expected;
    if (denom == 0.0) return 0.0;

    return (sum_nij - expected) / denom;
}
