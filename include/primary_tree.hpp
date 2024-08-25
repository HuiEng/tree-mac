
// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_primary_tree_HPP
#define INCLUDE_primary_tree_HPP
#include "tidy_tree.hpp"
// #include "temp_tree.hpp"

using namespace std;
size_t signatureWidth = signatureSize * sizeof(cell_type);

// typedef const cell_type *s_type;
typedef cell_type *s_type;
typedef const cell_type *const_s_type;
typedef vector<cell_type> sVec_type;

void printMatrtix(sVec_type temp_matrix)
{
    printMsg("printing matrices\n");
    for (size_t i = 0; i < temp_matrix.size(); i += signatureSize)
    {
        toBinaryIdx(stderr, &temp_matrix[i]);
    }
    printMsg("done print\n");
}

sVec_type createRandomSigs(const vector<cell_type> &sigs, size_t clusterCount, size_t s)
{
    unsigned seed = s;
    if (seed == 0)
    {
        seed = chrono::system_clock::now().time_since_epoch().count();
    }
    default_random_engine rng(seed);

    sVec_type clusterSigs(signatureSize * clusterCount);
    size_t signatureCount = sigs.size() / signatureSize;
    uniform_int_distribution<size_t> dist(0, signatureCount - 1);
    bool finished = false;

    unordered_set<string> uniqueSigs;
    for (size_t i = 0; i < signatureCount; i++)
    {
        size_t sig = dist(rng);
        string sigData(signatureWidth, ' ');
        memcpy(&sigData[0], &sigs[sig * signatureSize], signatureWidth);
        uniqueSigs.insert(sigData);
        if (uniqueSigs.size() >= clusterCount)
        {
            finished = true;
            break;
        }
    }

    size_t i = 0;
    for (const auto &sig : uniqueSigs)
    {
        memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureWidth);
        i++;
    }

    if (!finished)
    {
        if (uniqueSigs.size() != 1)
        {
            fprintf(stderr, "This should not happen\n");
            exit(1);
        }
        for (size_t i = 0; i < signatureSize; i++)
        {
            clusterSigs.push_back(clusterSigs[i]);
        }
    }

    return clusterSigs;
}

// set bit if at least have of the clusterSigs contains that bit
// finally, if result too little set bits, randomly set based on clusterSigs
sVec_type createMeanSig(const vector<cell_type> &clusterSigs)
{
    size_t seqCount = clusterSigs.size() / signatureSize;
    if (seqCount == 1)
    {
        return clusterSigs;
    }

    sVec_type meanSig(signatureSize);
    fill(&meanSig[0], &meanSig[0] + signatureSize, 0);
    if (seqCount == 0)
    {
        return meanSig;
    }

    vector<int> unflattenedSignature(signatureSize * bits_per_char);
    size_t c = 0;
    for (size_t i = 0; i < seqCount; i++)
    {
        const cell_type *signatureData = &clusterSigs[i * signatureSize];
        for (size_t i = 0; i < unflattenedSignature.size(); i++)
        {
            cell_type signatureMask = (cell_type)1 << (i % bits_per_char);

            if (signatureMask & signatureData[i / bits_per_char])
            {
                unflattenedSignature[i] += 1;
                c++;
            }
            // else
            // {
            //     unflattenedSignature[i] -= 1;
            // }
        }
    }
    size_t d = 0;

    cell_type *flattenedSignature = &meanSig[0];
    for (size_t i = 0; i < unflattenedSignature.size(); i++)
    {
        if (unflattenedSignature[i] >= (seqCount + 1) / 2)
        {
            d++;
            flattenedSignature[i / bits_per_char] |= (cell_type)1 << (i % bits_per_char);
        }
    }

    // size_t minSeqLength = countSetBits(&clusterSigs[0], signatureSize);
    // if (d < minSeqLength)
    // {
    //     // fprintf(stderr, "d: %zu, min %zu", d, minSeqLength);
    //     sVec_type randSig = createRandomSigs(clusterSigs, 1, 0);
    //     // fprintf(stderr, ", rand %zu", countSetBits(&randSig[0], signatureSize));
    //     for (size_t i = 0; i < signatureSize; i++)
    //     {
    //         meanSig[i] |= randSig[i];
    //     }
    //     // fprintf(stderr, ">>> %zu\n", countSetBits(&meanSig[0], signatureSize));

    // }

    return meanSig;
}

// Derived class
class primary_tree : public tidy_tree<s_type, const_s_type, sVec_type>
{
public:
    using tidy_tree::tidy_tree;

    inline void resizeMeans()
    {
        means.resize(capacity * signatureSize);
    }

    const_s_type getSeq(size_t i)
    {
        return &seqs[i * signatureSize];
    }

    const_s_type getSeq(sVec_type temp_seqs, size_t i)
    {
        return &temp_seqs[i * signatureSize];
    }

    s_type getMeanSig(size_t node) { return &means[node * signatureSize]; }

    // check if first child of branch is a singleton
    inline bool isSingleton(size_t child)
    {
        return matrices[child].size() == signatureSize;
    }

    void delSigFromMatrix(size_t node, size_t idx)
    {
        idx = idx * signatureSize;
        matrices[node].erase(matrices[node].begin() + idx, matrices[node].begin() + idx + signatureSize);
    }

    void readNodeSig(size_t parent, size_t child, const char *binFile)
    {
        vector<cell_type> signature;
        readSignatures(binFile, signature);

        updateMeanSig(child, &signature[0]);
        addSigToMatrix(parent, &signature[0]);
    }

    void readNodeMatrix(size_t child, const char *binFile)
    {
        readSignatures(binFile, matrices[child]);
    }

    void printSignature(ostream &wf, size_t node)
    {
        s_type mean_ = getMeanSig(node);
        wf.write((char *)mean_, signatureSize);
        // printMsg("node %zu\n", node);
        // toBinaryIdx(stderr, mean_);

        // // fprintf(stderr, ">>%zu\n", node);
        // // toBinaryIdx(stderr, mean_);
        // for (size_t i = 0; i < signatureSize; i++)
        // {
        //     wf.write(reinterpret_cast<const char *>(mean_ + i), sizeof(cell_type));
        // }
    }

    void printMatrix(ostream &wf, size_t node)
    {
        // s_type mean_ = &matrices[node][0];
        wf.write((char *)&matrices[node][0], signatureSize * seqIDs[node].size());
        // // fprintf(stderr, "node %zu, seq %zu\n", node, seqIDs[node].size());
        // for (size_t i =0; i<matrices[node].size();i+=signatureSize){
        //     s_type mean_ = &matrices[node][i];
        //     wf.write((char *)mean_, signatureSize);
        // }
    }

    void updateMeanSig(size_t node, const_s_type signature)
    {
        memcpy(&means[node * signatureSize], signature, signatureWidth);
    }

    ///?
    void addSigToMatrix(size_t node, const_s_type signature)
    {
        matrices[node].insert(matrices[node].end(), signature, signature + signatureSize);
    }

    double calcOverlapWrap(const_s_type a, const_s_type b)
    {
        return calcOverlap(a, b, signatureSize);
    }

    double calcSimilarityWrap(const_s_type a, const_s_type b)
    {
        return calcSimilarity(a, b, signatureSize);
    }

    size_t calcSetBitsWrap(const_s_type a)
    {
        return countSetBits(a, signatureSize);
    }

    //?
    double calcSimilaritySigToNode(size_t node, sVec_type signatures, size_t i)
    {
        return calcSimilarityWrap(&means[node * signatureSize], &signatures[i * signatureSize]);
    }

    sVec_type createNodeRandomSigs(size_t node, size_t clusterCount, size_t s)
    {
        return createRandomSigs(matrices[node], clusterCount, s);
    }
    

    inline sVec_type getNonAmbiMatrix(size_t node)
    {
        sVec_type temp_matrix;
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] != AMBI_T)
            {
                s_type sig = getMeanSig(child);
                temp_matrix.insert(temp_matrix.end(), sig, sig + signatureSize);
            }
        }

        return temp_matrix;
    }

    //?
    sVec_type unionNodeMean(size_t node)
    {
        sVec_type meanSig(signatureSize, 0);
        for (size_t child : childLinks[node])
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                meanSig[i] |= means[child * signatureSize + i];
            }
        }
        return meanSig;
    }

    sVec_type intersectNodeMean(size_t node)
    {
        sVec_type meanSig(signatureSize, 1);
        for (size_t child : childLinks[node])
        {
            for (size_t i = 0; i < signatureSize; i++)
            {
                meanSig[i] &= means[child * signatureSize + i];
            }
        }
        return meanSig;
    }

    inline int updateMatrixIdx(size_t node)
    {
        int idx = getNodeIdx(node);
        size_t parent = parentLinks[node];
        if (idx == -1)
        {
            fprintf(stderr, "ERROR: updateMatrixIdx cant find node %zu from parent %zu\n", node, parent);
            return idx;
        }
        memcpy(&matrices[parent][idx * signatureSize], getMeanSig(node), signatureWidth);
        return idx;
    }

    // // union mean of children
    // inline int updateNodeMean(size_t node)
    // {
    //     sVec_type meanSig;
    //     if (nodeType[node] == ROOT_T)
    //     {
    //         meanSig = unionNodeMean(node);
    //         // meanSig = createMeanSig(matrices[node]);
    //     }
    //     else if (nodeType[node] == SUPER_T)
    //     {
    //         meanSig = createMeanSig(getNonAmbiMatrix(node));
    //         // fprintf(stderr,"node: %zu, children: %zu, union bit: %zu, mean bit: %zu, inter bit: %zu\n",
    //         // node,
    //         // childLinks[node].size(),
    //         // calcSetBitsWrap(&unionNodeMean(node)[0]),
    //         // calcSetBitsWrap(&meanSig[0]),
    //         // calcSetBitsWrap(&intersectNodeMean(node)[0])
    //         // );
    //     }
    //     else if (nodeType[node] >= BRANCH_T)
    //     {
    //         meanSig = createMeanSig(getNonAmbiMatrix(node));
    //     }
    //     else
    //     {
    //         meanSig = createMeanSig(matrices[node]);
    //     }

    //     updateMeanSig(node, &meanSig[0]);
    //     updatePriority(node);
    //     return updateMatrixIdx(node);
    // }

    // union mean of children
    inline void updateNodeMean(size_t node)
    {
        sVec_type meanSig;
        if (nodeType[node] == ROOT_T)
        {
            meanSig = unionNodeMean(node);
            // meanSig = createMeanSig(matrices[node]);
        }
        else if (nodeType[node] == SUPER_T)
        {
            meanSig = createMeanSig(getNonAmbiMatrix(node));
            // fprintf(stderr,"node: %zu, children: %zu, union bit: %zu, mean bit: %zu, inter bit: %zu\n",
            // node,
            // childLinks[node].size(),
            // calcSetBitsWrap(&unionNodeMean(node)[0]),
            // calcSetBitsWrap(&meanSig[0]),
            // calcSetBitsWrap(&intersectNodeMean(node)[0])
            // );
        }
        else if (nodeType[node] >= BRANCH_T)
        {
            meanSig = createMeanSig(getNonAmbiMatrix(node));
        }
        else
        {
            meanSig = createMeanSig(matrices[node]);
        }

        updateMeanSig(node, &meanSig[0]);
        updatePriority(node);
    }

    inline void updateMatrixIdx(size_t parent, size_t idx, size_t node)
    {
        memcpy(&matrices[parent][idx * signatureSize], getMeanSig(node), signatureWidth);
    }

    // combine matrix of all descendant leaves of a node
    inline sVec_type createMeanSigAllMatrix(size_t node)
    {
        sVec_type tempMatrix;
        vector<size_t> leaves;
        getAllLeaves(node, leaves);
        for (size_t leaf : leaves)
        {
            insertVecRange(tempMatrix, matrices[leaf]);
        }
        sVec_type tempSig = createMeanSig(tempMatrix);

        // fprintf(stderr, "node: %zu, leaves: %zu,seqs: %zu, mean bit: %zu, temp bit: %zu, sim %.2f\n",
        //         node,
        //         leaves.size(),
        //         tempMatrix.size() / signatureSize,
        //         calcSetBitsWrap(getMeanSig(node)),
        //         calcSetBitsWrap(&tempSig[0]),
        //         calcSimilarityWrap(getMeanSig(node), &tempSig[0]));
        return tempSig;
    }

    // update nodes with all descendant matrix
    inline void updateAllMatrix(size_t node = 0)
    {
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] >= SUPER_T)
            {
                updateMeanSig(child, &createMeanSigAllMatrix(child)[0]);
                updateAllMatrix(child);
                // break;
            }
        }
    }

    // do after calling updateAllMatrix()
    // similar to prepReinsert except super is calculated based on all matrix
    inline void updateAllPriority(size_t node = 0)
    {
        for (size_t child : childLinks[node])
        {
            if (nodeType[child] >= SUPER_T)
            {
                sVec_type temp_matrix = matrices[child];
                size_t seqCount = temp_matrix.size() / signatureSize;
                double sumDistance = 0;
                for (size_t i = 0; i < temp_matrix.size(); i += signatureSize)
                {
                    double distance = calcSimilarityWrap(getMeanSig(child), &temp_matrix[i]); //?
                    sumDistance += distance;
                }
                double p = sumDistance / seqCount;
                // if (priority[child] != p)
                // {
                //     fprintf(stderr, "--- child %zu, before %.4f, after %.4f\n", child, priority[child], p);
                // }
                priority[child] = p;
                updateAllPriority(child);
            }
            else if (nodeType[child] <= LEAF_T)
            {
                updateNodeMean(child);
                seqIDs[child].clear();
                matrices[child].clear();
            }
            else
            { // branch
                updatePriority(child);
            }
        }
    }

    inline double calcNodeAvgSim(size_t node)
    {
        sVec_type temp_matrix = matrices[node];
        if (nodeType[node] >= BRANCH_T)
        {
            temp_matrix = getNonAmbiMatrix(node);
        }

        size_t seqCount = temp_matrix.size() / signatureSize;
        if (seqCount <= 1)
        {
            return 0;
        }

        double sumDistance = 0;
        s_type meanSig = getMeanSig(node);

        for (size_t i = 0; i < temp_matrix.size(); i += signatureSize)
        {
            double distance = calcSimilarityWrap(meanSig, &temp_matrix[i]); //?
            sumDistance += distance;
        }
        return sumDistance / seqCount;
    }

    void clearMean(size_t node)
    {
        s_type sig = getMeanSig(node);
        fill(sig, sig + signatureSize, 0);
    }

    void printNodeDistance(FILE *stream, sVec_type seqs, vector<size_t> &clusters)
    {
        fprintf(stream, "seq_id,clu,HD\n");
        for (size_t i = 0; i < clusters.size(); i++)
        {
            size_t tnode = clusters[i];
            double distance = calcHD(getMeanSig(tnode), &seqs[i * signatureSize]);
            fprintf(stream, "%zu,%zu,%.4f\n", i, tnode, distance);
        }
    }
};

#endif