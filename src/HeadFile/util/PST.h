#pragma once
#include <vector>
#include <algorithm>
#include <index/PPT.h>

struct PSTNode {
    PPT* ppt;
    int maxMu;
    PSTNode* left;
    PSTNode* right;
};

class PST2D {
private:
    PSTNode* root = nullptr;
    size_t nodeCount = 0; 
    PSTNode* build(std::vector<PPT*>& pts, int l, int r) {
        if (l > r) return nullptr;

        int mid = (l + r) >> 1;
        nodeCount++;
        PSTNode* node = new PSTNode{pts[mid], pts[mid]->mu, nullptr, nullptr};

        node->left  = build(pts, l, mid - 1);
        node->right = build(pts, mid + 1, r);

        if (node->left)  node->maxMu = std::max(node->maxMu, node->left->maxMu);
        if (node->right) node->maxMu = std::max(node->maxMu, node->right->maxMu);

        return node;
    }
    size_t memoryBytes() const {
        return nodeCount * sizeof(PSTNode);
    }
    void query(
        PSTNode* node,
        int epsQ, int muQ,
        std::vector<PPT*>& res
    ) const {
        if (!node || node->maxMu < muQ) return;

        if (node->ppt->epsilon >= epsQ && node->ppt->mu >= muQ)
            res.push_back(node->ppt);

        query(node->left, epsQ, muQ, res);
        query(node->right, epsQ, muQ, res);
    }

public:
    void build(std::vector<PPT>& ppts) {
        std::vector<PPT*> ptrs;
        for (auto& p : ppts) ptrs.push_back(&p);

        std::sort(ptrs.begin(), ptrs.end(),
            [](PPT* a, PPT* b) {
                return a->epsilon > b->epsilon;  // reverse
            });

        root = build(ptrs, 0, ptrs.size() - 1);
    }
    void printMemory() const {
        std::cout << "[PST] nodes = " << nodeCount
                  << ", memory = " << toMB(memoryBytes())
                  << " MB\n";
    }
    std::vector<PPT*> query(int epsQ, int muQ) const {
        std::vector<PPT*> res;
        query(root, epsQ, muQ, res);
        return res;
    }
};