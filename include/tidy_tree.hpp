#ifndef INCLUDE_TIDY_TREE_HPP
#define INCLUDE_TIDY_TREE_HPP
#include "similarity_status.hpp"

// Derived class
template <typename s_type, typename const_s_type, typename sVec_type>
class tidy_tree
{
public:
    size_t root = 0;                   // # of root node
    vector<int> nodeType;              // n entries, is this a branch node
    vector<size_t> potential_leaves;   // store seqIDs from singleton leaves
    vector<size_t> potential_sigs;     // n entries, is this a branch node
    vector<vector<size_t>> childLinks; // n * o entries, links to children
    vector<vector<size_t>> seqIDs;     // n * o entries, links to children
    vector<size_t> parentLinks;        // n entries, links to parents
    vector<distance_type> priority;    // n entries, links to parents
    // vector<vector<cell_type>> means;   // n * signatureSize entries, node signatures
    sVec_type means;            // n * signatureSize entries, node signatures
    vector<sVec_type> matrices; // capacity * signatureSize * n
    vector<omp_lock_t> locks;   // n locks
    size_t capacity = 0;        // Set during construction, currently can't change
    sVec_type seqs;

    virtual const_s_type getSeq(size_t i) { return returnEmpy<const_s_type>(); }
    virtual const_s_type getSeq(sVec_type temp_seqs, size_t i) { return returnEmpy<const_s_type>(); }

    virtual s_type getMeanSig(size_t node) { return returnEmpy<s_type>(); }

    virtual void updateAllMatrix(size_t node = 0) {}

    virtual void updateAllPriority(size_t node = 0) {}

    //*** need to do for sec tree
    void updateAll(size_t node = 0)
    {
        updateAllMatrix(node);
        updateAllPriority(node);
    }

    virtual inline bool isSingleton(size_t child) { return false; }

    virtual void testing(size_t node) {}

    virtual void delSigFromMatrix(size_t node, size_t idx) {}

    virtual void updateMeanSig(size_t node, const_s_type signature) {}

    virtual void addSigToMatrix(size_t node, const_s_type signature) {}

    virtual double calcSimilaritySigToNode(size_t node, sVec_type signatures, size_t i) { return 0; }

    virtual double calcSimilarityWrap(const_s_type a, const_s_type b) { return 0; }

    virtual size_t calcSetBitsWrap(const_s_type a) { return 0; }

    virtual double calcOverlapWrap(const_s_type a, const_s_type b) { return 0; }

    double calcScore(const_s_type a, const_s_type b, bool isRoot)
    {
        if (isRoot)
        {
            return calcOverlapWrap(a, b);
        }
        else
        {
            return calcSimilarityWrap(a, b);
        }
    }

    virtual sVec_type createNodeRandomSigs(size_t node, size_t clusterCount, size_t s = 0) { return returnEmpy<sVec_type>(); }

    virtual inline sVec_type getNonAmbiMatrix(size_t node) { return returnEmpy<sVec_type>(); }

    virtual sVec_type unionNodeMean(size_t node) { return returnEmpy<sVec_type>(); }

    // union mean of children
    virtual inline void updateNodeMean(size_t node) {}

    inline void updateLeafMatrix(size_t node)
    {
        matrices[node].clear();
        for (size_t idx : seqIDs[node])
        {
            addSigToMatrix(node, getSeq(idx));
        }
    }

    virtual inline void updateMatrixIdx(size_t parent, size_t idx, size_t node) {}

    virtual inline int updateMatrixIdx(size_t node) { return -1; }

    virtual inline double calcNodeAvgSim(size_t node) { return 0; }

    virtual void clearMean(size_t node) {}

    virtual inline void resizeMeans() {}

    void reserve(size_t capacity)
    {
        // For safety, only call this at startup currently
        if (this->capacity != 0)
        {
            fprintf(stderr, "Reserve can only be called from 0 capacity\n");
            exit(1);
        }
        this->capacity = capacity;

#pragma omp parallel
        {
#pragma omp single
            {
                nodeType.resize(capacity);
            }
#pragma omp single
            {
                childLinks.resize(capacity);
            }
#pragma omp single
            {
                seqIDs.resize(capacity);
            }
#pragma omp single
            {
                priority.resize(capacity);
            }
#pragma omp single
            {
                parentLinks.resize(capacity);
            }
#pragma omp single
            {
                locks.resize(capacity);
            }
#pragma omp single
            {
                matrices.resize(capacity);
            }
        }
    }

    tidy_tree(size_t capacity)
    {
        reserve(capacity);
        nodeType[root] = ROOT_T;
    }

    tidy_tree(size_t capacity, sVec_type input_seqs)
    {
        reserve(capacity);
        nodeType[root] = ROOT_T;
        seqs = input_seqs;
    }

    void printMatrix(FILE *stream, size_t node)
    {
        fprintf(stream, ">>>printing matrix at node %zu\n", node);
        for (s_type seq : matrices[node])
        {
            dbgPrintSignature(stream, seq);
        }
    }

    size_t getNewNodeIdx(vector<size_t> &insertionList)
    {
        if (insertionList.empty())
        {
            fprintf(stderr, "ERROR: ran out of insertion points\n");
            exit(1);
        }
        size_t idx = insertionList.back();
        insertionList.pop_back();

        // Initialise lock
        omp_init_lock(&locks[idx]);
        return idx;
    }

    size_t findAncestor(size_t node)
    {
        while (parentLinks[node] != root)
        {
            node = parentLinks[node];
        }
        return node;
    }

    tuple<size_t, size_t> findAncestorNlevel(size_t node)
    {
        size_t level = 0;
        while (parentLinks[node] != root)
        {
            node = parentLinks[node];
            level++;
        }
        return make_tuple(node, level);
    }

    size_t findLevel(size_t node)
    {
        size_t level = 0;
        while (node != root)
        {
            node = parentLinks[node];
            level++;
        }
        return level;
    }

    // get all descendant leaves of a node
    inline void getAllLeaves(size_t node, vector<size_t> &leaves)
    {
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] > LEAF_T)
            {
                getAllLeaves(child, leaves);
            }
            else if (nodeType[child] == LEAF_T)
            {
                leaves.push_back(child);
            }
        }
    }

    void printTreeJson_noRootLeaf(FILE *stream)
    {
        vector<size_t> candidates;
        for (size_t child : childLinks[root])
        {
            if (nodeType[child] > LEAF_T)
            {
                candidates.push_back(child);
            }
        }
        fprintf(stream, "var treeData = ");
        if (getChildCount(root) > 0)
        {
            printNodeJson(stream, root);
            fprintf(stream, "\",\"children\":[");

            printSubTreeJson(stream, candidates.back());
            candidates.pop_back();

            for (size_t child : candidates)
            {
                fprintf(stream, ",");
                printSubTreeJson(stream, child);
            }
            fprintf(stream, "]}");
        }
        fprintf(stream, ";\n");
    }

    void printNodeJson(FILE *stream, size_t tnode)
    {
        fprintf(stream, "{\"node\":\"%zu\",", tnode);
        fprintf(stream, "\"nodeType\":\"%d\",", nodeType[tnode]);
        fprintf(stream, "\"priority\":\"%.2f\",", priority[tnode]);
        fprintf(stream, "\"childCount\":\"%zu\",\"content\":\"*", seqIDs[tnode].size());
        for (size_t seq : seqIDs[tnode])
        {
            fprintf(stream, "%zu,", seq);
        }
        // fprintf(stream, "\",\"children\":[");
    }

    size_t getChildCount(size_t node)
    {
        return childLinks[node].size();
    }

    void printSubTreeJson(FILE *stream, size_t tnode)
    {
        if (getChildCount(tnode) > 0)
        {
            printNodeJson(stream, tnode);
            fprintf(stream, "\",\"children\":[");

            printSubTreeJson(stream, childLinks[tnode][0]);

            for (size_t i = 1; i < getChildCount(tnode); i++)
            {
                fprintf(stream, ",");
                printSubTreeJson(stream, childLinks[tnode][i]);
            }
            fprintf(stream, "]}");
        }
        else
        {
            // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, priority[tnode], seqIDs[tnode].size());
            // // fprintf(stream, "{\"node\":\"%zu\",\"priority\":\"%.2f\",\"childCount\":\"%zu\",\"content\":\"*", tnode, calcNodeMaxsimilarity(tnode), seqIDs[tnode].size());

            // for (size_t seq : seqIDs[tnode])
            // {
            //     fprintf(stream, "%zu,", seq);
            // }
            printNodeJson(stream, tnode);

            fprintf(stream, "\"}");
        }
    }

    void calcNodeDistance(FILE *stream, size_t tnode, size_t idx, const_s_type sig)
    {
        double distance = calcHD(getMeanSig(tnode), sig);
        fprintf(stream, "%zu,%zu,%.4f", tnode, idx, distance);
    }

    void printNodeDistance(FILE *stream, sVec_type seqs, vector<size_t> &clusters) {}

    void printTreeJson(FILE *stream, size_t node = 0)
    {
        fprintf(stream, "var treeData = ");
        printSubTreeJson(stream, node);
        fprintf(stream, ";\n");
    }

    virtual void printSignature(ostream &wf, size_t node) {}
    virtual void printMatrix(ostream &wf, size_t node) {}

    void printNodeMean(string outfolder, size_t node)
    {
        ofstream wf(outfolder + to_string(node) + ".bin", ios::out | ios::binary);
        writeInt(wf, signatureSize);
        // wf.write(reinterpret_cast<const char *>(&i), sizeof(cell_type));
        printSignature(wf, node);
        wf.close();
    }

    void printNodeMatrix(string outfolder, size_t node)
    {
        ofstream wf(outfolder + "m-" + to_string(node) + ".bin", ios::out | ios::binary);
        writeInt(wf, signatureSize);
        printMatrix(wf, node);
        wf.close();
    }

    void printSeqIDbyLeaf(FILE *stream, size_t leaf)
    {
        fprintf(stream, "%zu>", leaf);
        for (size_t id : seqIDs[leaf])
        {
            fprintf(stream, "%zu,", id);
        }
        fprintf(stream, "\n");
    }

    void printPotentialSigs(FILE *stream)
    {
        for (size_t id : potential_sigs)
        {
            fprintf(stream, "%zu,", id);
        }
        fprintf(stream, "\n");
    }

    void printLeavesWrap(FILE *stream, size_t node = 0)
    {

        for (size_t child : childLinks[node])
        {
            if (nodeType[child] <= LEAF_T)
            {
                printSeqIDbyLeaf(stream, child);
            }
            else
            {
                printLeavesWrap(stream, child);
            }
        }
    }

    void printLeaves(FILE *stream, size_t node = 0)
    {
        printPotentialSigs(stream);
        printLeavesWrap(stream, node);
    }

    void printTree(FILE *stream, string outfolder, size_t node)
    {
        fprintf(stream, "%zu", node);

        switch (nodeType[node])
        {
        case ROOT_T:
            fprintf(stream, "r");
            break;
        case SUPER_T:
            fprintf(stream, "s");
            break;
        case BRANCH_T:
            fprintf(stream, "b");
            break;
        case LEAF_T:
            fprintf(stream, "l");
            break;
        default:
            fprintf(stream, "a");
            break;
        }

        fprintf(stream, ">");
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] == AMBI_T)
            {
                fprintf(stream, "-");
            }
            fprintf(stream, "%zu,", child);
        }
        fprintf(stream, "\n");
        for (size_t child : childLinks[node])
        {
            if (getChildCount(child) > 0)
            {
                printTree(stream, outfolder, child);
            }
        }
    }

    void printTree(string outfolder, vector<size_t> &insertionList)
    {
        FILE *tFile = fopen((outfolder + "tree.txt").c_str(), "w");
        fprintf(tFile, "%zu\n%zu\n", signatureSize, getNewNodeIdx(insertionList));
        printTree(tFile, outfolder, 0);
    }

    virtual void readNodeSig(size_t parent, size_t child, const char *binFile) {}
    virtual void readNodeMatrix(size_t child, const char *binFile) {}

    void readNode(size_t parent, size_t child)
    {
        parentLinks[child] = parent;
        childLinks[parent].push_back(child);
    }

    int getNodeIdx(size_t node)
    {
        size_t parent = parentLinks[node];
        int idx = -1;
        for (size_t i = 0; i < getChildCount(parent); i++)
        {
            if (childLinks[parent][i] == node)
            {
                idx = i;
                break;
            }
        }
        if (idx == -1)
        {
            fprintf(stderr, "ERROR finding %zu from %zu!!\n", node, parent);
        }
        return idx;
    }

    s_type createRandomSig(vector<size_t> children, size_t s = 0)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        unsigned seed = s;
        if (seed == 0)
        {
            seed = chrono::system_clock::now().time_since_epoch().count();
        }
        default_random_engine rng(seed);
        uniform_int_distribution<size_t> dist(0, children.size() - 1);

        const_s_type randomSig; // = getMeanSig(children[dist(rng)])
        // find the smallest windows count
        size_t winNum = getMeanSig(children[0]).size();
        for (size_t child : children)
        {
            s_type matrix = getMeanSig(child);
            if (matrix.size() < winNum)
            {
                winNum = matrix.size();
            }
        }

        randomSig.resize(winNum);

        vector<vector<size_t>> counters;
        for (size_t w = 0; w < winNum; w++)
        {
            size_t s = dist(rng);
            randomSig[w] = getMeanSig(children[s])[w];
        }

        return randomSig;
    }

    // delete a child from its parent, need to format the child separately
    inline void deleteNode(size_t node)
    {
        // printMsg("deleting %zu\n", node);
        size_t parent = parentLinks[node];

        int idx = getNodeIdx(node);

        if (idx != -1)
        {
            // removeVecIdx(matrices[parent], idx);
            delSigFromMatrix(parent, idx);
            removeVecIdx(childLinks[parent], idx);
        }
    }

    void clearNode(size_t node)
    {
        nodeType[node] = 0;
        childLinks[node].clear();
        seqIDs[node].clear();
        matrices[node].clear();
        clearMean(node);
        priority[node] = 0;
        parentLinks[node] = 0;
    }

    size_t deleteUnitig(size_t node)
    {
        int idx = getNodeIdx(node);
        if (idx < 0)
        {
            return root;
        }
        size_t parent = parentLinks[node];
        size_t child = childLinks[node][0];

        childLinks[parent][idx] = child;
        parentLinks[child] = parent;
        // deleteNode(node);
        clearNode(node);

        return parent;
    }

    inline bool checkdeleteUniSuper(size_t node)
    {
        if (nodeType[node] != ROOT_T && getChildCount(node) == 1)
        {
            // deleteUnitig(node);
            // return true;
            return (deleteUnitig(node) == root);
        }
        return false;
    }

    inline void updateParentMean(size_t node)
    {
        while (node != root)
        {
            updateNodeMean(node);
            if (updateMatrixIdx(node) == -1)
            {
                break;
            }
            node = parentLinks[node];
        }
    }

    // for updateTree, only only parent of the last child
    inline void updateParentMean_last(size_t node)
    {
        while (node != root)
        {
            updateNodeMean(node);
            size_t parent = parentLinks[node];
            if (updateMatrixIdx(node) != getChildCount(parent) - 1)
            {
                break;
            }
            node = parent;
        }
    }

    // union mean of children
    inline void updatePriority(size_t node)
    {
        priority[node] = calcNodeAvgSim(node);
    }

    inline size_t stayNode(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        //? if (isBranchNode[node])
        if (nodeType[node] > LEAF_T)
        {
            return tt_root(signature, insertionList, idx, node);
        }
        else
        {
            seqIDs[node].push_back(idx);
            addSigToMatrix(node, signature);
            // updatePriority(node);
            return node;
        }
    }

    void moveParent(size_t child, size_t new_parent, bool d = true)
    {
        childLinks[new_parent].push_back(child);
        addSigToMatrix(new_parent, getMeanSig(child));
        // printMsg("moving %zu to %zu \n", child, new_parent);

        if (d)
        {
            deleteNode(child);
        }
        parentLinks[child] = new_parent;
        // printMsg("deleted %zu \n", child);
    }

    // haven't update childLinks[t_parent] yet
    inline size_t createParent(size_t node, vector<size_t> &insertionList)
    {
        size_t t_parent = getNewNodeIdx(insertionList);
        nodeType[t_parent] = BRANCH_T;
        parentLinks[t_parent] = node;

        // update node
        childLinks[node].push_back(t_parent);

        printMsg("create t_parent %zu\n", t_parent);
        return t_parent;
    }

    inline size_t createNode(const_s_type signature, vector<size_t> &insertionList, size_t t_parent, size_t idx)
    {
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        addSigToMatrix(new_node, signature);
        seqIDs[new_node].push_back(idx);
        updateMeanSig(new_node, signature);
        nodeType[new_node] = LEAF_T;

        // then add new node to t_parent
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, getMeanSig(new_node));

        if (nodeType[t_parent] >= BRANCH_T)
        {
            printMsg("--------------b %zu, %.2f\n", new_node, t_parent);
            updatePriority(t_parent);
        }

        printMsg("create new node %zu\n", new_node);
        return new_node;
    }

    // Ambi Node cannot be the first child or else update node mean will have problem
    inline size_t createAmbiNode(const_s_type signature, vector<size_t> &insertionList, size_t node, size_t idx)
    {
        size_t dest = createNode(signature, insertionList, node, idx);

        printMsg(">>Ambi\n");
        nodeType[dest] = AMBI_T;
        updatePriority(node);
        return dest;
    }

    inline pair<size_t, double> similarityStatusF(const_s_type sig1, const_s_type sig2, double local_stay_t, double offset = 1, bool isBranch = false)
    {
        //?
        if (local_stay_t > stay_threshold)
        {
            local_stay_t = stay_threshold;
        }

        double local_split_t = split_threshold * offset;
        double similarity = calcSimilarityWrap(sig1, sig2);

        printMsg("%.2f (%.2f)", similarity, local_split_t);

        if (similarity >= local_stay_t)
        {
            printMsg(": STAY>\n");
            // return STAY_F;
            return make_pair(STAY_F, similarity);
        }
        else if (similarity <= local_split_t)
        {
            printMsg(": SPLIT>\n");
            // return SPLIT_F;
            return make_pair(SPLIT_F, similarity);
        }
        else
        {
            printMsg(": NN>\n");
            // return NN_F;
            if (isBranch)
            {
                return make_pair(NN_BRANCH_F, similarity);
            }
            else
            {
                return make_pair(NN_LEAVE_F, similarity);
            }
        }
    }

    inline pair<size_t, double> similarityStatus(size_t child, const_s_type signature)
    {
        printMsg("<%zu, ", child);

        if (nodeType[child] == SUPER_T)
        {
            // super can only be stay or split
            printMsg("(super %.2f)", priority[child]);
            // return similarityStatusF(getMeanSig(child), signature, priority[child], 100);
            return similarityStatusF(getMeanSig(child), signature, split_threshold, 100);
            // return similarityStatusF(getMeanSig(child), signature, split_node_threshold, 100);
        }
        else if (nodeType[child] == BRANCH_T)
        {
            // super can only be stay or split

            if (priority[child] == 0)
            {
                // branch only 1 child
                printMsg("(branch)");
                return similarityStatusF(getMeanSig(child), signature, stay_threshold, 1, true);
            }
            else
            {
                printMsg("(branch %.2f)", priority[child]);
                return similarityStatusF(getMeanSig(child), signature, priority[child], 1, true);
            }
        }
        else
        {
            return similarityStatusF(getMeanSig(child), signature, stay_threshold);
        }
    }

    inline void mergeBranch(size_t branch, vector<size_t> &insertionList)
    {
        size_t parent = parentLinks[branch];
        vector<size_t> children = childLinks[branch];
        clearNode(branch);
        parentLinks[branch] = parent;
        nodeType[branch] = LEAF_T;

        printMsg("merging branch %zu\n", branch);
        for (size_t child : children)
        {
            if (nodeType[child] != AMBI_T)
            {
                insertVecRange(matrices[branch], matrices[child]);
                insertVecRange(seqIDs[branch], seqIDs[child]);
            }
            clearNode(child);
        }

        if (matrices[branch].size() == 0)
        {
            fprintf(stderr, "ERROR merging branch\n");
            removeAmbi(insertionList, branch);
        }
        else
        {
            updateParentMean(branch);
        }
    }

    inline vector<size_t> getAmbiChildren(size_t node)
    {
        vector<size_t> ambiChildren;
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] == AMBI_T)
            {
                ambiChildren.push_back(child);
            }
        }
        return ambiChildren;
    }

    inline size_t insertAmbi(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t branch)
    {
        vector<size_t> ambiChildren = getAmbiChildren(branch);
        if (ambiChildren.size() == 0)
        {
            return createAmbiNode(signature, insertionList, branch, idx);
        }

        printMsg("InsertAmbi %zu\n", branch);

        size_t dest = 0;
        double max_similarity = 0;
        for (size_t ambi : ambiChildren)
        {

            double similarity = calcSimilarityWrap(getMeanSig(ambi), signature);
            printMsg(">>ambi %zu, %.2f\n", ambi, similarity);
            if (similarity >= stay_threshold)
            {
                dest = stayNode(signature, insertionList, idx, ambi);
                nodeType[dest] = LEAF_T;
                updateParentMean(dest);
                printMsg("@@ambi %zu, %.2f\n", dest, priority[ambi]);
                return dest;
            }
            else if (similarity > max_similarity)
            {
                max_similarity = similarity;
                dest = ambi;
            }
        }

        if (dest == 0)
        {
            printMsg("??? ambi %zu\n", branch);
            return createAmbiNode(signature, insertionList, branch, idx);
        }
        else
        {
            return stayNode(signature, insertionList, idx, dest);
        }
    }

    // create ambi and insert to it if there is none
    // else, add to suitable ambi and turn the ambi to leaf if the prioroty >stay threshold
    // else create new ambi
    inline size_t insertBranch(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t branch)
    {
        double max_similarity = stay_threshold;
        size_t dest = 0;
        for (size_t child : childLinks[branch])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg("b %zu, %.2f\n", child, similarity);
            if (similarity >= max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }

        if (dest != 0)
        {
            dest = stayNode(signature, insertionList, idx, dest);
            if (nodeType[dest] == AMBI_T)
            {
                nodeType[dest] = LEAF_T;
                updateParentMean(dest);
            }
            return dest;
        }

        return insertAmbi(signature, insertionList, idx, branch);
    }

    // find best leaf from the NN branches, if dest = 0, create new node
    inline size_t doNNBranch(const_s_type signature, size_t node, tt_data dt, vector<size_t> &insertionList, size_t idx)
    {
        // xxx, should be ok
        size_t t_parent = node;
        vector<size_t> target = dt.nn_branch;
        if (checkCreateSuper(node, true, dt.mismatch, target.size()))
        {
            t_parent = createSuper(node, insertionList, target);
            updateParentMean(t_parent);
        }

        printMsg("doNNBranch %zu\n", t_parent);

        double max_similarity = stay_threshold;
        size_t dest = 0;
        double max_ambi_similarity = split_threshold;
        size_t dest_ambi = 0;
        for (size_t branch : target)
        {
            for (size_t child : childLinks[branch])
            {
                double similarity = calcSimilarityWrap(getMeanSig(child), signature);
                if (nodeType[child] == AMBI_T)
                {
                    printMsg("b_ambi %zu, %.2f\n", child, similarity);
                    if (similarity >= max_ambi_similarity)
                    {
                        max_ambi_similarity = similarity;
                        dest_ambi = child;
                    }
                }
                else
                {
                    printMsg("b %zu, %.2f\n", child, similarity);
                    if (similarity >= max_similarity)
                    {
                        max_similarity = similarity;
                        dest = child;
                    }
                }
            }
        }

        if (dest + dest_ambi == 0)
        {
            dest = createNode(signature, insertionList, t_parent, idx);
        }
        else if (dest != 0)
        {
            dest = stayNode(signature, insertionList, idx, dest);
        }
        else
        {
            dest = stayNode(signature, insertionList, idx, dest_ambi);
            updateNodeMean(dest_ambi);
            if (priority[dest_ambi] >= stay_threshold)
            {
                printMsg("@@ambi %zu, %.2f\n", dest_ambi, priority[dest_ambi]);
                nodeType[dest_ambi] = LEAF_T;
                // updateParentMean(dest_ambi);// will do below
            }
        }

        updateParentMean(dest);
        return dest;
    }

    inline size_t createUniBranch(size_t node, vector<size_t> &insertionList, const_s_type signature, size_t idx)
    {
        size_t t_parent = createParent(node, insertionList);
        size_t dest = createNode(signature, insertionList, t_parent, idx);
        addSigToMatrix(node, signature);
        updateParentMean(t_parent);
        return dest;
    }

    // create a new branch under "node"
    // then move candidates to the new branch
    inline size_t createBranch(size_t node, vector<size_t> &insertionList, vector<size_t> candidates)
    {
        size_t t_parent = createParent(node, insertionList);

        for (size_t child : candidates)
        {
            moveParent(child, t_parent);
        }

        addSigToMatrix(node, getMeanSig(t_parent));
        updateParentMean(t_parent);
        return t_parent;
    }

    // create a new super under "node"
    // then move candidates to the new super
    inline size_t createSuper(size_t node, vector<size_t> &insertionList, vector<size_t> candidates)
    {
        size_t t_parent = createBranch(node, insertionList, candidates);
        nodeType[t_parent] = SUPER_T;
        return t_parent;
    }

    //? return if all children moved
    inline bool relocate(vector<size_t> parents, vector<size_t> children)
    {
        size_t moved = 0;
        printMsg("Relocating...\n");
        for (size_t child : children)
        {
            double max_simimlarity = split_threshold;
            size_t dest = 0;
            for (size_t parent : parents)
            {
                double similarity = calcSimilarityWrap(getMeanSig(child), getMeanSig(parent));
                printMsg("reloc %zu, %.2f\n", parent, similarity);
                if (similarity > max_simimlarity)
                {
                    max_simimlarity = similarity;
                    dest = parent;
                }
            }
            //? maybe find the best level to move
            if (dest != 0 && dest != parentLinks[child])
            {
                moveParent(child, dest);
                moved++;
            }
        }

        if (moved > 0)
        {
            for (size_t node : parents)
            {
                updateNodeMean(node);
            }
        }
        return moved == children.size();
    }

    // separate roots and non roots children of node
    // first list contains roots
    // second list contains non-roots
    inline vector<vector<size_t>> separateRootChildren(size_t node)
    {
        vector<vector<size_t>> children(2);
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] == ROOT_T)
            {
                children[0].push_back(child);
            }
            else
            {
                children[1].push_back(child);
            }
        }
        return children;
    }

    inline tt_data getStatus(const_s_type signature, size_t node)
    {
        tt_data dt;

        for (size_t child : childLinks[node])
        {
            // skip ambiNode
            if (nodeType[child] == AMBI_T)
            {
                // printMsg("Ambi: %zu\n", child);
                continue;
            }

            double similarity = 0;
            size_t status = 0;

            if (nodeType[child] == ROOT_T)
            { ///???
                double overlap = calcOverlapWrap(getMeanSig(child), signature);
                similarity = calcSimilarityWrap(getMeanSig(child), signature);
                printMsg("(root %zu, %.2f, %.2f)\n", child, priority[child], similarity);
                if (overlap > stay_threshold)
                {
                    if (similarity > dt.max_root_similarity)
                    {
                        dt.max_root_similarity = similarity;
                        dt.dest_root = child;
                    }
                    dt.statuses.set(7);
                }
                continue;
            }
            std::tie(status, similarity) = similarityStatus(child, signature);
            switch (status)
            {
            case STAY_F:
                // stay.push_back(child);
                if (nodeType[child] == SUPER_T)
                {
                    dt.stay_super.push_back(child);
                    if (similarity > dt.max_super_similarity)
                    {
                        dt.max_super_similarity = similarity;
                        dt.dest_super = child;
                    }
                    dt.statuses.set(5);
                }
                else if (nodeType[child] == BRANCH_T)
                {
                    dt.stay_branch.push_back(child);
                    if (similarity > dt.max_branch_similarity)
                    {
                        dt.max_branch_similarity = similarity;
                        dt.dest_branch = child;
                    }
                    dt.statuses.set(3);
                }
                else
                {
                    dt.stay_leaf.push_back(child);
                    if (similarity > dt.max_similarity)
                    {
                        dt.max_similarity = similarity;
                        dt.dest = child;
                    }
                    // printMsg("child %zu\n", dt.dest);
                    dt.statuses.set(1);
                }

                break;
            case SPLIT_F:
                // mismatch.push_back(child);
                dt.mismatch++;
                break;
            case NN_LEAVE_F:
                dt.nn_leaf.push_back(child);
                dt.statuses.set(2);
                break;
            default:
                dt.nn_branch.push_back(child);
                dt.statuses.set(4);
                break;
            }
        }
        return dt;
    }

    // return true if it will not create a unitigrelocateRemaining(node,
    inline bool checkCreateSuper(size_t node, bool moved_all, size_t mismatch, size_t target_size)
    {
        bool doRoot = false;
        bool doSuper = false;

        if (nodeType[node] == ROOT_T)
        {
            if (target_size > 1)
            {
                doRoot = true;
            }
            else if (!moved_all)
            {
                doRoot = true;
            }
        }
        else
        {
            if (!moved_all && mismatch != 0)
            {
                doSuper = true;
            }
            else if (target_size > 1 && mismatch != 0)
            {
                doSuper = true;
            }
        }

        return doRoot || doSuper;
    }

    inline bool relocateRemaining(size_t node, bool moved_all_current, vector<size_t> candidates, vector<size_t> targets)
    {
        bool moved_all_candidates = moved_all_current;
        if (!moved_all_candidates)
        {
            // try to relocate the remaining stay leaves to super
            vector<size_t> remaining_candidates;
            for (size_t child : candidates)
            {
                if (parentLinks[child] == node)
                {
                    remaining_candidates.push_back(child);
                }
            }

            moved_all_candidates = relocate(targets, remaining_candidates);
        }
        return moved_all_candidates;
    }

    inline void movedRemaining(size_t node, size_t t_parent, bool moved_all, vector<size_t> candidates)
    {
        if (!moved_all)
        {
            for (size_t child : candidates)
            {
                if (parentLinks[child] == node)
                {
                    moveParent(child, t_parent);
                }
            }
        }
    }

    // relocate candidates to target, then createSuper accordingly
    // level = 3 = super [default]
    // level = 2 = nn branch
    // level = 1 = stay branch
    inline size_t doTarget_candidates(tt_data dt, size_t node, vector<size_t> candidates, vector<size_t> &insertionList, size_t level = 3)
    {
        vector<size_t> target;
        switch (level)
        {
        case 1:
            target = dt.stay_branch;
            break;
        case 2:
            target = dt.nn_branch;
            break;
        case 3:
            target = dt.stay_super;
            break;
        default:
            fprintf(stderr, "ERROR - doTarget_StayNN\n");
            return node;
            break;
        }
        bool moved_all_nn = relocate(target, candidates);

        if (!checkdeleteUniSuper(node))
        {
            if (checkCreateSuper(node, moved_all_nn, dt.mismatch, target.size()))
            { // there is some nn branches can't be moved to stay super, create new super
                size_t t_parent = createSuper(node, insertionList, target);
                movedRemaining(node, t_parent, moved_all_nn, candidates);
                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    // relocate stay_candidates and nn_candidates to target, then createSuper accordingly
    // level = 3 = super [default]
    // level = 2 = nn branch
    // level = 1 = stay branch
    inline size_t doTarget_StayNN(tt_data dt, size_t node, vector<size_t> stay_candidates, vector<size_t> nn_candidates, vector<size_t> &insertionList, size_t level = 3)
    {
        vector<size_t> target;
        switch (level)
        {
        case 1:
            target = dt.stay_branch;
            break;
        case 2:
            target = dt.nn_branch;
            break;
        case 3:
            target = dt.stay_super;
            break;
        default:
            fprintf(stderr, "ERROR - doTarget_StayNN\n");
            return node;
            break;
        }

        bool moved_all_stay = relocate(target, stay_candidates);
        bool moved_all_nn = relocate(target, nn_candidates);
        bool moved_all = moved_all_stay && moved_all_nn;

        if (checkdeleteUniSuper(node))
        {
            return node;
        }

        if (checkCreateSuper(node, moved_all, dt.mismatch, target.size()))
        { // there is some stay and/or nn candidates can't be moved to stay X, create new Super
            size_t t_parent = createSuper(node, insertionList, target);
            movedRemaining(node, t_parent, moved_all_stay, stay_candidates);
            movedRemaining(node, t_parent, moved_all_nn, nn_candidates);
            updateParentMean(t_parent);
            return t_parent;
        }

        return node;
    }

    // move low to high
    // then remaining low to target (stay_super)
    // then move high to target
    // step 2&3 reversible, same effect
    inline size_t doStaySuperLowHigh(tt_data dt, size_t node, vector<size_t> low_candidates, vector<size_t> high_candidates, vector<size_t> &insertionList)
    {
        vector<size_t> target = dt.stay_super;
        bool moved_all_low = relocate(high_candidates, low_candidates);
        moved_all_low = relocateRemaining(node, moved_all_low, low_candidates, target);
        bool moved_all_high = relocate(target, high_candidates);
        bool moved_all = moved_all_low && moved_all_high;

        if (!checkdeleteUniSuper(node))
        {
            if (checkCreateSuper(node, moved_all, dt.mismatch, target.size()))
            { // there is some nn highes can't be moved to stay super, create new super
                size_t t_parent = createSuper(node, insertionList, target);
                movedRemaining(node, t_parent, moved_all_low, low_candidates);
                movedRemaining(node, t_parent, moved_all_high, high_candidates);
                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    // relocate stay (level 1, default) or nn leaf (level 2) to stay branch, or else stay super, or else nn branch
    // then stay branch to super
    // then nn branch to super
    inline size_t doLeaf_others(tt_data dt, size_t node, vector<size_t> &insertionList, size_t level = 1)
    {
        vector<size_t> candidates;
        if (level == 1)
        {
            candidates = dt.stay_leaf;
        }
        else
        {
            candidates = dt.nn_leaf;
        }

        bool moved_all_leaf = relocate(dt.stay_branch, candidates);
        moved_all_leaf = relocateRemaining(node, moved_all_leaf, candidates, dt.stay_super);
        moved_all_leaf = relocateRemaining(node, moved_all_leaf, candidates, dt.nn_branch);
        bool moved_all_stay_branch = relocate(dt.stay_super, dt.stay_branch);
        bool moved_all_nn_branch = relocate(dt.stay_super, dt.nn_branch);
        bool moved_all = moved_all_leaf && moved_all_stay_branch && moved_all_nn_branch;

        if (!checkdeleteUniSuper(node))
        {
            if (checkCreateSuper(node, moved_all, dt.mismatch, candidates.size()))
            { // there is some nn highes can't be moved to stay super, create new super
                size_t t_parent = createSuper(node, insertionList, candidates);
                movedRemaining(node, t_parent, moved_all_leaf, candidates);
                movedRemaining(node, t_parent, moved_all_stay_branch, dt.stay_branch);
                movedRemaining(node, t_parent, moved_all_nn_branch, dt.nn_branch);

                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    // relocate candidates to high, then remaining to low
    inline size_t doCandidates_HighLow(tt_data dt, size_t node, vector<size_t> candidates, vector<size_t> high_parents, vector<size_t> low_parents, vector<size_t> &insertionList)
    {
        bool moved_all = relocate(high_parents, candidates);
        moved_all = relocateRemaining(node, moved_all, candidates, low_parents);

        if (!checkdeleteUniSuper(node))
        {
            if (nodeType[node] == ROOT_T || dt.mismatch != 0)
            {
                // create super for high and low
                size_t t_parent = createSuper(node, insertionList, high_parents);
                for (size_t child : low_parents)
                {
                    moveParent(child, t_parent);
                }

                // there is some candidates can't be moved to high or low, move to super
                movedRemaining(node, t_parent, moved_all, candidates);
                updateParentMean(t_parent);
                return t_parent;
            }
        }
        return node;
    }

    inline size_t doBit1(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;
        size_t pos;
        for (pos = 0; pos < dt.statuses.size(); ++pos)
        {
            if (dt.statuses.test(pos))
            {
                break;
            }
        }

        switch (pos)
        {
        case 1:
            printMsg(">>Stat>> Stay Leaf-\n");
            if (dt.stay_leaf.size() > 1)
            {
                createBranch(node, insertionList, dt.stay_leaf);
            }
            return stayNode(signature, insertionList, idx, dt.dest);
            break;
        case 2:
            printMsg(">>Stat>> NN Leaf-\n");
            t_branch = createBranch(node, insertionList, dt.nn_leaf);
            return insertBranch(signature, insertionList, idx, t_branch);
            // return createAmbiNode(signature, insertionList, t_branch, idx);
            break;
        case 3:
            printMsg(">>Stat>> Stay Branch-\n");
            if (dt.stay_branch.size() > 1)
            {
                if (nodeType[node] == ROOT_T || dt.mismatch > 0)
                {
                    t_parent = createSuper(node, insertionList, dt.stay_branch);
                }
            }
            return insertBranch(signature, insertionList, idx, dt.dest_branch);
            break;
        case 4:
            printMsg(">>Stat>> NN Branch- %zu, mismatch:%zu\n", dt.nn_branch.size(), dt.mismatch);

            if (dt.nn_branch.size() > 1)
            {

                return doNNBranch(signature, node, dt, insertionList, idx);
            }
            else
            {
                return insertBranch(signature, insertionList, idx, dt.nn_branch[0]);
            }
            break;
        case 5:
            printMsg(">>Stat>> Stay Super-\n");
            return tt_super(signature, insertionList, idx, dt.dest_super);
            break;
        case 7:
            printMsg(">>Stat>> Stay Root-\n");
            return tt_root(signature, insertionList, idx, dt.dest_root);
            break;
        default:
            fprintf(stderr, "ERROR 1 bit idx %zu, pos %zu!!\n", idx, pos);
            return 0;
            break;
        }
    }

    inline size_t doBit2(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;
        size_t tempbits = (dt.statuses & bitset<8>("00101010")).count();
        // only stay
        if (tempbits == 2)
        {
            if (dt.statuses.test(5))
            {
                if (dt.statuses.test(1))
                {
                    printMsg(">>Stat>> Stay Leaf and Stay Super-\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    doTarget_candidates(dt, node, dt.stay_leaf, insertionList);
                }
                else
                {
                    printMsg(">>Stat>> Stay Branch and Stay Super-\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    doTarget_candidates(dt, node, dt.stay_branch, insertionList);
                }
            }
            else
            {
                printMsg(">>Stat>> Stay Leaf and Stay Branch-\n");
                dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                doTarget_candidates(dt, node, dt.stay_leaf, insertionList, 1);
            }
            return dt.dest;
        }
        else if (tempbits == 1)
        {
            if (dt.statuses.test(4))
            {
                if (dt.statuses.test(1))
                {
                    printMsg(">>Stat>> Stay Leaf and NN Branch-\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    doTarget_candidates(dt, node, dt.stay_leaf, insertionList, 2);
                    return dt.dest;
                }
                else if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> Stay Branch and NN Branch-\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    if (nodeType[node] == ROOT_T || dt.mismatch != 0)
                    {
                        t_parent = createSuper(node, insertionList, dt.stay_branch);
                        for (size_t b : dt.nn_branch)
                        {
                            moveParent(b, t_parent);
                        }
                        updateParentMean(t_parent);
                    }
                    return dt.dest;
                }
                else if (dt.statuses.test(5))
                { //?
                    printMsg(">>Stat>> NN Branch and Stay Super-\n");
                    // doTarget_candidates(dt, node, dt.nn_branch, insertionList);
                    return tt_super(signature, insertionList, idx, dt.dest_super);
                }
            }
            else
            {
                // printMsg("NN Leaf and ");
                if (dt.statuses.test(1))
                {
                    printMsg(">>Stat>> Stay Leaf and NN Leaf-\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    t_branch = createBranch(node, insertionList, dt.stay_leaf);
                    for (size_t l : dt.nn_leaf)
                    {
                        moveParent(l, t_branch);
                    }
                    return dt.dest;
                }
                else if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> NN Leaf and Stay Branch-\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    doTarget_candidates(dt, node, dt.nn_leaf, insertionList, 1);
                    return dt.dest;
                }
                else if (dt.statuses.test(5))
                {
                    printMsg(">>Stat>> NN Leaf and Stay Super-\n");
                    // doTarget_candidates(dt, node, dt.nn_leaf, insertionList);
                    return tt_super(signature, insertionList, idx, dt.dest_super);
                }
            }
        }
        else
        {
            printMsg(">>Stat>> NN Leaf and NN Branch-\n");
            t_parent = node;
            if (dt.mismatch != 0)
            {
                t_parent = createSuper(node, insertionList, dt.nn_branch);
                for (size_t leaf : dt.nn_leaf)
                {
                    moveParent(leaf, t_parent);
                }
            }

            dt.dest = createNode(signature, insertionList, t_parent, idx);
            updateParentMean(t_parent);
            return dt.dest;
        }

        fprintf(stderr, "ERROR 2 bit %zu!!", idx);
        return 0;
    }

    inline size_t doBit3(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;
        if (dt.statuses.test(5))
        {
            if (dt.statuses.test(1))
            {
                dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                if (dt.statuses.test(2))
                {
                    printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Super-\n");
                    doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_leaf, insertionList);
                }
                else if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> Stay Leaf and Stay Branch and Stay Super-\n");
                    doStaySuperLowHigh(dt, node, dt.stay_leaf, dt.stay_branch, insertionList);
                }
                else
                {
                    printMsg(">>Stat>> Stay Leaf and NN branch and Stay Super-\n");
                    doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_branch, insertionList);
                }

                return dt.dest;
            }
            else if (dt.statuses.test(2))
            {
                // printMsg("Stay Super and NN Leaf and ");
                if (dt.statuses.test(3))
                {
                    printMsg(">>Stat>> NN Leaf and Stay Branch and Stay Super-\n");
                    dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                    doStaySuperLowHigh(dt, node, dt.nn_leaf, dt.stay_branch, insertionList);
                    return dt.dest;
                }
                else
                {
                    printMsg(">>Stat>> NN Leaf and NN branch and Stay Super-\n");
                    // doStaySuperLowHigh(dt, node, dt.nn_leaf, dt.nn_branch, insertionList);
                    return tt_super(signature, insertionList, idx, dt.dest_super);
                }
            }
            else
            {
                printMsg(">>Stat>> Stay Branch and NN branch and Stay Super-\n");
                dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                doTarget_StayNN(dt, node, dt.stay_branch, dt.nn_branch, insertionList);
                return dt.dest;
            }
        }
        else
        {
            if (dt.statuses.test(1))
            {
                if (dt.statuses.test(2))
                {
                    if (dt.statuses.test(3))
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch-\n");
                        dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                        doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_leaf, insertionList, 1);
                    }
                    else
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and NN Branch-\n");
                        dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                        doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_leaf, insertionList, 2);
                    }
                }
                else
                {
                    printMsg(">>Stat>> Stay Leaf and Stay Branch and NN branch-\n");
                    dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                    doCandidates_HighLow(dt, node, dt.stay_leaf, dt.stay_branch, dt.nn_branch, insertionList);
                }
            }
            else
            {
                printMsg(">>Stat>> NN Leaf and Stay Branch and NN branch-\n");
                dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
                doCandidates_HighLow(dt, node, dt.nn_leaf, dt.stay_branch, dt.nn_branch, insertionList);
            }

            return dt.dest;
        }

        fprintf(stderr, "ERROR 3 bit %zu!!", idx);
        return 0;
    }

    inline size_t doBit4(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t t_parent, t_branch;

        if (dt.statuses.test(1))
        {
            dt.dest = stayNode(signature, insertionList, idx, dt.dest);
            if (dt.statuses.test(2))
            {
                insertVecRange(dt.stay_leaf, dt.nn_leaf);
                if (dt.statuses.test(3))
                {
                    if (dt.statuses.test(4))
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch-\n");
                        doCandidates_HighLow(dt, node, dt.stay_leaf, dt.stay_branch, dt.nn_branch, insertionList);
                    }
                    else
                    {
                        printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch and Stay Super-\n");
                        doStaySuperLowHigh(dt, node, dt.stay_leaf, dt.stay_branch, insertionList);
                    }
                }
                else
                {
                    printMsg(">>Stat>> Stay Leaf and NN Leaf and NN Branch and Stay Super-\n");
                    doTarget_StayNN(dt, node, dt.stay_leaf, dt.nn_branch, insertionList);
                }
            }
            else
            {
                printMsg(">>Stat>> Stay Leaf and Stay Branch and NN Branch and Stay Super-\n");
                dt.dest = stayNode(signature, insertionList, idx, dt.dest);
                doLeaf_others(dt, node, insertionList);
            }
            return dt.dest;
        }
        else
        {
            printMsg(">>Stat>> NN Leaf and Stay Branch and NN Branch and Stay Super-\n");
            dt.dest = insertBranch(signature, insertionList, idx, dt.dest_branch);
            doLeaf_others(dt, node, insertionList, 2);
            return dt.dest;
        }

        fprintf(stderr, "ERROR 4 bit %zu!!", idx);
        return 0;
    }

    inline tt_data getStatus_potential(const_s_type signature)
    {
        tt_data dt;
        for (size_t i : potential_sigs)
        {
            double similarity = calcSimilarityWrap(getSeq(i), signature);

            if (similarity >= stay_threshold)
            {
                dt.stay_leaf.push_back(i);
                dt.statuses.set(1);
            }
            else if (similarity <= split_threshold)
            {
                dt.mismatch++;
            }
            else
            {
                dt.nn_leaf.push_back(i);
                dt.statuses.set(2);
            }
        }
        return dt;
    }

    // dest cannot be super or root
    // in the case of stay size = 1
    // dest = stay[0]
    inline size_t growtree_without_root(tt_data dt, const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        if (dt.statuses.none())
        {
            // reserve for insert_potential
            if (nodeType[node] == ROOT_T)
            {
                return node;
            }
            return createNode(signature, insertionList, node, idx);
        }

        size_t t_parent, t_branch;
        size_t setbits = dt.statuses.count();
        if (setbits == 1)
        {
            return doBit1(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 2)
        {
            return doBit2(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 3)
        {
            return doBit3(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 4)
        {
            return doBit4(dt, signature, insertionList, idx, node);
        }
        else if (setbits == 5)
        {
            printMsg(">>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch and Stay Super-\n");
            t_branch = createBranch(node, insertionList, dt.stay_leaf);
            for (size_t l : dt.nn_leaf)
            {
                moveParent(l, t_branch);
            }
            printMsg("combined leaves to branch %zu\n", t_branch);
            dt.stay_branch.push_back(t_branch);
            dt.dest_branch = t_branch;
            dt.statuses.reset(1);
            dt.statuses.reset(2);
            return doBit3(dt, signature, insertionList, idx, node);
        }

        printMsg(">>Stat>> ERROR!!\n");
        return 0;
    }

    inline size_t tt(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {
        tt_data dt = getStatus(signature, node);
        return growtree_without_root(dt, signature, insertionList, idx, node);
    }

    inline size_t getStatus(size_t idx)
    {
        const_s_type signature = getSeq(idx);
        return getStatusFlag(getStatus(signature, 0));
    }

    inline size_t tt_super(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node)
    {
        size_t dest = searchBestNonRoot(signature, node);
        return stayNode(signature, insertionList, idx, dest);
    }

    inline void dissolveSuper(size_t parent)
    {
        printMsg(">>> Dissolving super %zu\n", parent);
        size_t grandparent = parentLinks[parent];
        vector<size_t> candidates = childLinks[parent];
        for (size_t child : candidates)
        {
            moveParent(child, grandparent);
        }
        deleteNode(parent);
        clearNode(parent);
        updateParentMean(grandparent);
    }

    inline size_t findBestRoot(const_s_type signature, size_t node)
    {
        double max_similarity = 0;
        double max_overlap = 0;
        size_t best_root = node;
        for (size_t child : childLinks[node])
        {
            double overlap = calcOverlapWrap(getMeanSig(child), signature);
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg("(root %zu, %.2f, %.2f)\n", child, overlap, similarity);
            if (overlap > max_overlap)
            {
                max_overlap = overlap;
                max_similarity = similarity;
                best_root = child;
            }
            else if (overlap == max_overlap)
            {
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    best_root = child;
                }
            }
        }

        if (best_root == node)
        {
            best_root = childLinks[node][0];
            fprintf(stderr, "ERROR: Search Best Root - force search node %zu => child %zu\n", node, best_root);
        }
        return best_root;
    }

    inline size_t tt_root(const_s_type signature, vector<size_t> &insertionList, size_t idx, size_t node = 0)
    {

        if (getChildCount(node) == 0 || nodeType[childLinks[node][0]] != ROOT_T)
        {
            // omp_set_lock(&locks[node]);
            size_t dest = tt(signature, insertionList, idx, node);
            // omp_set_lock(&locks[node]);
            // updateParentMean(dest);
            return dest;
        }

        size_t best_root = findBestRoot(signature, node);
        size_t dest = 0;
        printMsg("traverse root %zu\n", best_root);
        dest = tt_root(signature, insertionList, idx, best_root);

        if (dest == 0)
        {
            fprintf(stderr, "tt_root %zu\n", idx);
        }

        addSubtree(insertionList, best_root);
        return dest;
    }

    inline size_t first_insert(vector<size_t> &insertionList, size_t idx, size_t parent = 0)
    {
        return createNode(getSeq(idx), insertionList, root, idx);
    }

    inline size_t insert(vector<size_t> &insertionList, size_t idx)
    {
        const_s_type signature = getSeq(idx);
        size_t node = tt(signature, insertionList, idx);
        if (node == 0)
        {
            node = insert_potential(insertionList, idx);
        }
        printMsg("inserted %zu at %zu\n\n", idx, node);
        return node;
    }

    inline size_t insert_potential(vector<size_t> &insertionList, size_t idx, size_t t_parent = 0)
    {
        const_s_type signature = getSeq(idx);
        tt_data dt = getStatus_potential(signature);
        if (dt.statuses.none())
        {
            potential_sigs.push_back(idx);
            return 0;
        }

        size_t dest = createNode(signature, insertionList, t_parent, idx);
        printMsg("***>>Stat>> Stay Leaf- %zu\n", dt.stay_leaf.size(), dt.nn_leaf.size());
        // combine stay and NN for potential
        insertVecRange(dt.stay_leaf, dt.nn_leaf);
        for (size_t i : dt.stay_leaf)
        {

            addSigToMatrix(dest, getSeq(i));
            removeVecValue(potential_sigs, i);
        }
        insertVecRange(seqIDs[dest], dt.stay_leaf);
        return dest;

        // size_t setbits = dt.statuses.count();
        // size_t before = potential_sigs.size();

        // if (setbits == 1)
        // {
        //     if (dt.statuses.test(1))
        //     {
        //         size_t dest = createNode(signature, insertionList, t_parent, idx);
        //         printMsg("***>>Stat>> Stay Leaf- %zu\n", dt.stay_leaf.size());
        //         for (size_t i : dt.stay_leaf)
        //         {

        //             addSigToMatrix(dest, getSeq(i));
        //             removeVecValue(potential_sigs, i);
        //         }
        //         insertVecRange(seqIDs[dest], dt.stay_leaf);
        //         return dest;
        //     }
        //     else if (dt.statuses.test(2))
        //     {
        //         printMsg("***>>Stat>> NN Leaf-%zu\n", dt.nn_leaf.size());
        //         size_t branch = createParent(t_parent, insertionList);
        //         size_t dest = createNode(signature, insertionList, branch, idx);
        //         size_t i = dt.nn_leaf.back();
        //         dt.nn_leaf.pop_back();

        //         size_t ambiNode = createAmbiNode(getSeq(i), insertionList, branch, i);
        //         removeVecValue(potential_sigs, i);
        //         for (size_t i : dt.nn_leaf)
        //         {
        //             addSigToMatrix(ambiNode, getSeq(i));
        //             removeVecValue(potential_sigs, i);
        //         }
        //         insertVecRange(seqIDs[ambiNode], dt.nn_leaf);
        //         updateNodeMean(branch);
        //         addSigToMatrix(t_parent, getMeanSig(branch));
        //         return dest;
        //     }
        // }
        // else if (dt.statuses.test(1) && dt.statuses.test(2))
        // {
        //     printMsg("***>>Stat>> Stay Leaf and NN Leaf-%zu,%zu\n", dt.stay_leaf.size(), dt.nn_leaf.size());
        //     size_t branch = createParent(t_parent, insertionList);
        //     size_t dest = createNode(signature, insertionList, branch, idx);
        //     for (size_t i : dt.stay_leaf)
        //     {
        //         addSigToMatrix(dest, getSeq(i));
        //         removeVecValue(potential_sigs, i);
        //     }
        //     insertVecRange(seqIDs[t_parent], dt.stay_leaf);

        //     size_t i = dt.nn_leaf.back();
        //     dt.nn_leaf.pop_back();

        //     size_t ambiNode = createAmbiNode(getSeq(i), insertionList, branch, i);
        //     removeVecValue(potential_sigs, i);
        //     for (size_t i : dt.nn_leaf)
        //     {
        //         addSigToMatrix(ambiNode, getSeq(i));
        //         removeVecValue(potential_sigs, i);
        //     }
        //     insertVecRange(seqIDs[ambiNode], dt.nn_leaf);
        //     updateNodeMean(branch);
        //     addSigToMatrix(t_parent, getMeanSig(branch));
        //     return dest;
        // }

        // fprintf(stderr, "ERROR: this should not happen in insert_potential: seq %zu\n", idx);
        // return 0;
    }

    size_t kMeans(size_t node, vector<vector<size_t>> &clusters, size_t clusterCount)
    {
        sVec_type temp_centroids = createNodeRandomSigs(node, clusterCount, kMean_seed);

        for (size_t child : childLinks[node])
        {
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < clusterCount; i++)
            {
                double similarity = calcSimilaritySigToNode(child, temp_centroids, i);
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            clusters[dest].push_back(child);
        }

        for (vector<size_t> cluster : clusters)
        {
            if (cluster.size() <= 1)
            {
                return 0;
            }
        }

        return 1;
    }

    size_t kMeansLeaf(size_t node, vector<vector<size_t>> &clusters, size_t clusterCount)
    {
        sVec_type temp_centroids = createNodeRandomSigs(node, clusterCount, kMean_seed);

        for (size_t idx : seqIDs[node])
        {
            size_t dest = 0;
            double max_similarity = 0;
            for (size_t i = 0; i < clusterCount; i++)
            {
                double similarity = calcSimilarityWrap(getSeq(temp_centroids, i), getSeq(idx));
                if (similarity > max_similarity)
                {
                    max_similarity = similarity;
                    dest = i;
                }
            }
            clusters[dest].push_back(idx);
        }

        for (vector<size_t> cluster : clusters)
        {
            if (cluster.size() <= 1)
            {
                return 0;
            }
        }

        return 1;
    }

    size_t addSubtree(vector<size_t> &insertionList, size_t node, size_t clusterCount = 2)
    {
        if (childLinks[node].size() <= tree_order)
        {
            return 0;
        }

        vector<vector<size_t>> clusters(clusterCount);
        if (kMeans(node, clusters, clusterCount) == 0)
        {
            return 0;
        }

        printMsg(">> addSubtree %zu\n", node);

        size_t grandparent = parentLinks[node];
        // omp_set_lock(&locks[grandparent]);

        size_t t_parent = createParent(grandparent, insertionList);
        nodeType[t_parent] = ROOT_T;
        for (size_t child : clusters[1])
        {
            moveParent(child, t_parent);
        }
        updateNodeMean(t_parent);
        addSigToMatrix(grandparent, getMeanSig(t_parent));

        updateParentMean(node);
        // omp_unset_lock(&locks[grandparent]);

        return 1;
    }

    size_t splitLeaf(size_t node, vector<size_t> &insertionList)
    {
        if (nodeType[node] != LEAF_T)
        {
            return 0;
        }

        // if (seqIDs[node].size() < 10)
        // {
        //     return 0;
        // }

        if (priority[node] > split_threshold)
        {
            return 0;
        }

        size_t clusterCount = 2;
        vector<vector<size_t>> clusters(clusterCount);
        if (kMeansLeaf(node, clusters, clusterCount) == 0)
        {
            fprintf(stderr, ">> tried splitLeaf %zu\n", node);
            return 0;
        }

        // do original leaf
        seqIDs[node] = clusters[0];
        updateLeafMatrix(node);

        // split to new sibling
        size_t t_parent = parentLinks[node];
        size_t new_node = getNewNodeIdx(insertionList);
        parentLinks[new_node] = t_parent;
        seqIDs[new_node] = clusters[1];
        updateLeafMatrix(new_node);
        nodeType[new_node] = LEAF_T;
        updateNodeMean(new_node);
        childLinks[t_parent].push_back(new_node);
        addSigToMatrix(t_parent, getMeanSig(new_node));

        fprintf(stderr, ">> splitLeaf %zu to new node %zu\n", node, new_node);
        return 1;
    }

    size_t forceSplitRoot(vector<size_t> &insertionList, size_t node = 0, size_t clusterCount = 2)
    {
        if (getChildCount(node) <= tree_order)
        {
            return 0;
        }
        vector<vector<size_t>> clusters(clusterCount);
        if (kMeans(node, clusters, clusterCount) == 0)
        {
            return 0;
        }

        printMsg(">> forceSplitRoot %zu\n", node);

        // reuse clusterSize to store the new t_parents
        for (vector<size_t> cluster : clusters)
        {
            size_t t_parent = createParent(node, insertionList);
            nodeType[t_parent] = ROOT_T;
            for (size_t child : cluster)
            {
                moveParent(child, t_parent);
            }
            updateNodeMean(t_parent);
            addSigToMatrix(node, getMeanSig(t_parent));
        }
        updateParentMean(node);

        return 1;
    }

    inline size_t insertSplitRoot(vector<size_t> &insertionList, size_t idx)
    {
        const_s_type signature = getSeq(idx);
        size_t node = tt_root(signature, insertionList, idx);
        if (nodeType[node] == ROOT_T)
        {
            node = insert_potential(insertionList, idx, node);
        }
        printMsg("inserted %zu at %zu,%zu\n\n", idx, node);
        forceSplitRoot(insertionList);
        return node;
    }

    inline size_t searchBestNonRoot(const_s_type signature, size_t node = 0)
    {
        printMsg("\nSearch Best Non Root %zu\n ", node);
        size_t dest = 0;
        double max_similarity = 0;
        // double avg_similarity = 0;
        vector<size_t> candidates;

        for (size_t child : childLinks[node])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg(" (%zu, %.2f) ", child, similarity);

            if (similarity >= max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }
        printMsg("\n>>>%zu\n", dest);

        if (nodeType[dest] > LEAF_T)
        {
            return searchBestNonRoot(signature, dest);
        }
        else
        {
            return dest;
        }
    }

    inline size_t searchAllNonRoot(const_s_type signature, size_t node = 0)
    {
        printMsg("\nSearch All Non Root %zu\n ", node);
        size_t dest = 0;
        double max_similarity = 0;

        for (size_t child : childLinks[node])
        {

            if (nodeType[child] > LEAF_T)
            {
                double similarity = calcSimilarityWrap(getMeanSig(child), signature);
                if (similarity <= split_threshold)
                {
                    continue;
                }
                child = searchAllNonRoot(signature, child);
            }
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg(" (%zu, %.2f) ", child, similarity);

            if (similarity >= max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }
        printMsg("\n>>>%zu\n", dest);

        return dest;
    }

    // grandchildren may be mix, use searchBest() to check
    inline size_t searchBestRoot(const_s_type signature, size_t node = 0)
    {
        printMsg("\nSearch Best Root %zu\n ", node);

        if (getChildCount(node) == 0)
        {
            // this shouldn't happen;
            fprintf(stderr, "ERROR: Search Best Root for root %zu has no child\n", node);
        }

        double max_similarity = 0;
        double max_overlap = 0;
        size_t dest = findBestRoot(signature, node);

        printMsg("\n>>>%zu\n", dest);

        if (nodeType[dest] > LEAF_T)
        {
            return searchBest(signature, dest); // should always be this
        }
        else
        {
            return dest;
        }
    }

    // decide if children contains roots only, non-roots only or mix
    // for mix, solve non roots first then solve roots
    // compare result for both and return best
    inline size_t searchBest(const_s_type signature, size_t node = 0)
    {
        vector<vector<size_t>> children = separateRootChildren(node);
        // if no root
        if (children[0].size() == 0)
        {
            return searchBestNonRoot(signature, node);
        }
        else if (children[1].size() == 0)
        {
            return searchBestRoot(signature, node);
        }
        printMsg("\nSearch Best Mix %zu\n ", node);

        double dest_similarity = 0;

        // do non root, include super, branch and leaf
        size_t dest = 0;
        double max_similarity = 0;
        for (size_t child : children[1])
        {
            double similarity = calcSimilarityWrap(getMeanSig(child), signature);
            printMsg(" (%zu, %.2f) ", child, similarity);

            if (similarity >= max_similarity)
            {
                max_similarity = similarity;
                dest = child;
            }
        }

        if (nodeType[dest] > LEAF_T)
        {
            dest = searchBestNonRoot(signature, dest);
        }

        dest_similarity = calcSimilarityWrap(getMeanSig(dest), signature);

        // early pruning
        if (dest_similarity >= stay_threshold)
        {
            return dest;
        }

        // do root
        size_t best_root = 0;
        double max_similarity_root = 0;
        double max_overlap_root = 0;

        for (size_t root_ : children[0])
        {
            double overlap = calcOverlapWrap(getMeanSig(root_), signature);
            double similarity = calcSimilarityWrap(getMeanSig(root_), signature);
            printMsg(" (r%zu, %.2f, %.2f) ", root_, overlap, similarity);

            if (overlap > max_overlap_root)
            {
                max_overlap_root = overlap;
                max_similarity_root = similarity;
                best_root = root_;
            }
            else if (overlap == max_overlap_root)
            {
                if (similarity > max_similarity_root)
                {
                    max_similarity_root = similarity;
                    best_root = root_;
                }
            }
        }

        if (best_root == 0)
        {
            // fprintf(stderr, "ERROR: search best mix\n");
            return dest;
        }

        if (nodeType[best_root] > LEAF_T)
        {
            best_root = searchBest(signature, best_root);
        }

        printMsg("\n>>>best_root: %zu, dest: %zu\n", best_root, dest);
        max_overlap_root = calcOverlapWrap(getMeanSig(best_root), signature);
        if (max_overlap_root > dest_similarity)
        {
            return best_root;
        }
        else
        {
            return dest;
        }
    }

    // cluster means still remain
    void prepReinsert(size_t node = 0)
    {
        //?
        if (nodeType[node] <= LEAF_T)
        {
            updateNodeMean(node);
            seqIDs[node].clear();
            matrices[node].clear();
        }
        else
        {
            updatePriority(node);
        }
        // fprintf(stderr, "<node %zu, type %zu, set %zu, priority %.4f>\n", node, nodeType[node], calcSetBitsWrap(getMeanSig(node)), priority[node]);

        for (size_t child : childLinks[node])
        {
            prepReinsert(child);
        }
    }

    void updateTree(size_t node = 0)
    {
        if (nodeType[node] == AMBI_T)
        {
            // updatePriority(node);
        }
        if (nodeType[node] == LEAF_T)
        {
            updateParentMean_last(node);
        }
        else
        {
            for (size_t child : childLinks[node])
            {
                updateTree(child);
            }
        }
    }

    size_t reinsert(size_t idx)
    {
        const_s_type signature = getSeq(idx);
        size_t node = searchBest(signature);
        double dest_similarity = calcSimilarityWrap(getMeanSig(node), signature);
        omp_set_lock(&locks[node]);
        // if (dest_similarity <= 0.1)
        // {
        //     fprintf(stderr, "%zu,%zu,%.4f\n", idx, node, dest_similarity);
        // }
        seqIDs[node].push_back(idx);
        addSigToMatrix(node, signature);
        omp_unset_lock(&locks[node]);
        return node;
    }

    inline size_t deleteParentUnitig(size_t parent)
    {
        while (getChildCount(parent) <= 1 && parent != root)
        {
            // fprintf(stderr, "---deleteUnitig %zu\n", parent);
            if (getChildCount(parent) == 0)
            {
                size_t node = parent;
                parent = parentLinks[node];
                deleteNode(node);
                clearNode(node);
            }
            else
            {
                parent = deleteUnitig(parent);
            }
        }

        return parent;
    }

    // // get all descendant leaves with seqCount less than {singleton} and ambi of a node
    // void removeRedundant(size_t node = 0)
    // {
    //     vector<size_t> children = childLinks[node];
    //     for (size_t child : children)
    //     {
    //         if (nodeType[child] == AMBI_T)
    //         {
    //             // fprintf(stderr, "*delete %zu,%zu,%zu\n", child, nodeType[child], seqIDs[child].size());
    //             deleteNode(child);
    //             clearNode(child);
    //         }
    //         else if (nodeType[child] == LEAF_T)
    //         {
    //             if (seqIDs[child].size() <= singleton)
    //             {
    //                 // fprintf(stderr, "*delete %zu,%zu,%zu\n", child, nodeType[child], seqIDs[child].size());
    //                 deleteNode(child);
    //                 clearNode(child);
    //             }
    //             else
    //             {
    //                 updateNodeMean(child);
    //             }
    //         }
    //         else if (nodeType[child] == BRANCH_T && priority[child] >= stay_threshold)
    //         {
    //             mergeBranch(child);
    //         }
    //         else
    //         {
    //             removeRedundant(child);
    //         }
    //     }
    //     size_t parent = deleteParentUnitig(node);
    //     if (parent != root)
    //     {
    //         updateNodeMean(parent);
    //     }
    // }

    void removeAmbi(vector<size_t> &insertionList, size_t node = 0)
    {
        size_t cleared = 0;

        vector<size_t> children = childLinks[node];
        if (nodeType[node] >= SUPER_T)
        {
            for (size_t child : children)
            {
                removeAmbi(insertionList, child);
            }
        }
        else if (nodeType[node] == BRANCH_T)
        {
            if (priority[node] >= stay_threshold)
            {
                mergeBranch(node, insertionList);
            }
            else
            {
                for (size_t child : children)
                {
                    removeAmbi(insertionList, child);
                }
            }
        }
        else if (nodeType[node] == AMBI_T)
        {
            // fprintf(stderr, "*delete %zu,%zu,%zu\n", node, nodeType[node], seqIDs[node].size());
            deleteNode(node);
            cleared = node;
        }
        else if (seqIDs[node].size() <= singleton)
        {
            insertVecRange(potential_leaves, seqIDs[node]);
            deleteNode(node);
            cleared = node;
        }
        else
        {
            updateParentMean(node);
            // splitLeaf(node, insertionList);
        }
        deleteParentUnitig(parentLinks[node]);
        if (cleared != 0)
        {
            clearNode(cleared);
        }
    }

    void doSmallLeave(vector<size_t> &insertionList)
    {
        fprintf(stderr, "doSmallLeave %zu\n", potential_leaves.size());
        for (size_t i : potential_leaves)
        {
            insert(insertionList, i);
        }
        potential_leaves.clear();
    }

    void doSmallLeave_split(vector<size_t> &insertionList)
    {
        fprintf(stderr, "doSmallLeave_split %zu\n", potential_leaves.size());
        for (size_t i : potential_leaves)
        {
            insertSplitRoot(insertionList, i);
        }
        potential_leaves.clear();
    }

    void outputHierarchy(FILE *pFile, size_t node = 0, size_t rank = 0)
    {
        for (size_t child : childLinks[node])
        {
            fprintf(pFile, "%zu,%zu,%zu\n", node, child, rank);
            if (nodeType[child] > LEAF_T)
            {
                outputHierarchy(pFile, child, rank + 1);
            }
        }
    }

    void destroyLocks(size_t node)
    {
        omp_destroy_lock(&locks[node]);
        if (nodeType[node] > LEAF_T)
        {
            for (size_t child : childLinks[node])
            {
                destroyLocks(child);
            }
        }
    }

    inline void destroyLocks()
    {
        destroyLocks(root);
    }
};

#endif