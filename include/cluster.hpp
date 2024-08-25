#ifndef INCLUDE_CLUSTER_HPP
#define INCLUDE_CLUSTER_HPP
#include <omp.h>
using namespace std;

bool random_ = false;
size_t cap = 0;
size_t iteration = 0;
bool iteration_given = false;
bool debug_ = false;
bool force_split_ = false;
bool skip_ = false;
unsigned seed;

struct tree_meta_st
{
    bool readTree_ = false;
    bool writeTree_ = false;
    const char *inputTreePath;
    const char *outputTreePath;
    bool writeJSON = false;
};

tree_meta_st tree_meta;

template <typename cmdline_type>
string setArgs(cmdline_type args)
{
    split_threshold = 0.5;
    stay_threshold = 0.8;
    minimiser_match_threshold = 3;

    //?
    if (!args.input_given)
    {
        cout << "No input and/or query given! Exiting...\n";
        exit(1);
    }

    string inputFile = args.input_arg;

    // if (args.skip_arg)
    // {
    //     skip_ = true;
    //     fprintf(stderr, "Skip inital tree construction, thresholds have no use\n");
    // }

    if (args.singleton_given)
    {
        singleton = args.singleton_arg;
        fprintf(stderr, "Min cluster size = %zu\n", singleton);
    }

    if (args.topology_in_given && !args.topology_out_given)
    {
        skip_ = true;
        fprintf(stderr, "Skip inital tree construction, thresholds have no use\n");
    }

    if (args.random_arg)
    {
        random_ = true;
        seed = args.seed_arg;
    }

    if (args.iteration_given)
    {
        iteration_given = true;
        iteration = args.iteration_arg;
    }

    if (args.stay_threshold_given)
    {
        stay_threshold = args.stay_threshold_arg;
    }

    if (args.split_threshold_given)
    {
        split_threshold = args.split_threshold_arg;
    }
    else if (args.topology_in_given)
    {
        // continue;
    }
    else if (args.single_arg)
    {
        if (args.multiple_arg)
        {
            split_threshold = getSplitThresholdListSingle(inputFile);
        }
        else
        {
            split_threshold = getSplitThresholdSingle(inputFile);
        }
    }
    else
    {
        split_threshold = getSplitThreshold(inputFile, args.multiple_arg);
    }

    split_node_threshold = split_threshold / 2;

    fprintf(stderr, "split threshold: %.2f\n", split_threshold);
    fprintf(stderr, "stay threshold: %.2f\n", stay_threshold);
    fprintf(stderr, "split_node_threshold threshold: %.2f\n", split_node_threshold);

    if (args.minimiser_match_given)
    {
        minimiser_match_threshold = args.minimiser_match_arg;
        fprintf(stderr, "minimiser_match threshold: %zu\n", minimiser_match_threshold);
    }

    cap = args.sizeCap_arg;

    if (args.capacity_given)
    {
        partree_capacity = args.capacity_arg;
        fprintf(stderr, "partree_capacity: %zu\n", partree_capacity);
    }

    debug_ = args.debug_arg;
    print_ = args.print_arg;
    updateAll_ = args.updateAll_arg;
    force_split_ = args.force_split_arg;
    kMean_seed = args.kMean_seed_arg;
    fprintf(stderr, "kMean_seed: %zu\n", kMean_seed);
    if (args.tree_order_given)
    {
        tree_order = args.tree_order_arg;
        fprintf(stderr, "tree_order: %zu\n", tree_order);
    }

    tree_meta.readTree_ = args.topology_in_given;
    tree_meta.writeTree_ = args.topology_out_given;
    tree_meta.inputTreePath = args.topology_in_arg;
    tree_meta.outputTreePath = args.topology_out_arg;
    tree_meta.writeJSON = args.print_tree_given;

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    auto fileName = rawname + "-s" + to_string((int)(stay_threshold * 100)) + "-l" + to_string((int)(split_threshold * 100));
    if (args.tag_given)
    {
        fileName = args.tag_arg + fileName;
    }

    return (fileName + ".txt");
}

template <typename tree_type>
size_t loadSeqIDs(const string folder, tree_type &tree)
{
    string s;
    size_t leaves = 0;
    size_t leaf = 0;
    string delimiter = ",";

    // Read from the text file
    ifstream listStream((folder + "clusters.txt").c_str());

    while (getline(listStream, s))
    {

        size_t pos = s.find(">");
        string node = s.substr(0, pos);
        sscanf(node.c_str(), "%zu", &leaf);
        s.erase(0, pos + 1);
        string seqIDStr = "";
        size_t seqID = 0;
        set<size_t> parents;
        while ((pos = s.find(delimiter)) != string::npos)
        {
            seqIDStr = s.substr(0, pos);
            sscanf(seqIDStr.c_str(), "%zu", &seqID);
            s.erase(0, pos + 1);

            tree.seqIDs[leaf].push_back(seqID);
        }

        size_t parent = tree.parentLinks[leaf];
        parents.insert(parent);

        if (tree.nodeType[leaf] != AMBI_T)
        {
            tree.nodeType[leaf] = LEAF_T;
            // tree.readNodeMatrix(leaf, (folder + "m-" + to_string(leaf) + ".bin").c_str());
            for (size_t i : tree.seqIDs[leaf])
            {
                tree.addSigToMatrix(leaf, tree.getSeq(i));
            }
            tree.updateNodeMean(leaf);
        }
        tree.addSigToMatrix(parent, tree.getSeq(tree.seqIDs[leaf][0]));

        do
        {
            set<size_t> parents_temp = parents;
            parents_temp.erase(0);
            parents.clear();
            for (size_t parent : parents_temp)
            {
                tree.updateNodeMean(parent);
                size_t grandparent = tree.parentLinks[parent];
                tree.addSigToMatrix(grandparent, tree.getMeanSig(parent));
                parents.insert(grandparent);
            }

        } while (parents.size() != 0);
        leaves++;
    }

    // Close the file
    listStream.close();

    return leaves;
}

template <typename tree_type>
string readTreeLine(string s, string folder, tree_type &tree)
{
    char node_type = 0;
    size_t parent = 0;
    size_t child = 0;
    string delimiter = ">";

    size_t pos = s.find(delimiter);
    string node = s.substr(0, pos);

    sscanf(node.c_str(), "%zu%s", &parent, &node_type);
    switch (node_type)
    {
    case 'r':
        tree.nodeType[parent] = ROOT_T;
        break;
    case 's':
        tree.nodeType[parent] = SUPER_T;
        break;
    case 'b':
        tree.nodeType[parent] = BRANCH_T;
        break;
    default:
        break;
    }
    s.erase(0, pos + 1);

    delimiter = ",";
    string childStr = "";
    while ((pos = s.find(delimiter)) != string::npos)
    {
        bool isAmbi = false;
        childStr = s.substr(0, pos);
        if (childStr[0] == '-')
        {
            isAmbi = true;
            childStr = s.substr(1, pos - 1);
        }
        sscanf(childStr.c_str(), "%zu", &child);
        s.erase(0, pos + 1);

        tree.readNode(parent, child);
        if (isAmbi)
        {
            tree.nodeType[child] = AMBI_T;
        }
    }
    return node;
}

// output signature size, if file too big => return empty "sigs", use readSignaturesBatch
template <typename tree_type>
size_t readTree(const string folder, tree_type &tree)
{
    string line;

    // Read from the text file
    ifstream listStream((folder + "tree.txt").c_str());

    // read signatureSize
    getline(listStream, line);

    // read last nodeIdx
    getline(listStream, line);
    size_t offset = 0;
    sscanf(line.c_str(), "%zu", &offset);
    if (offset > partree_capacity)
    {
        fprintf(stderr, "ERROR: partree_capacity must be greater than %zu\n", offset + 1);
        exit(1);
    }

    while (getline(listStream, line))
    {
        readTreeLine(line, folder, tree);
    }

    // Close the file
    listStream.close();

    fprintf(stderr, "Loaded %zu nodes\n", offset);
    return offset;
}

void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
    }
}

void compressClusterList(vector<size_t> &clusters)
{
    unordered_map<size_t, size_t> remap;
    for (size_t &clus : clusters)
    {
        if (remap.count(clus))
        {
            clus = remap[clus];
        }
        else
        {
            size_t newClus = remap.size();
            remap[clus] = newClus;
            clus = newClus;
        }
    }
    fprintf(stderr, "Output %zu clusters\n", remap.size());
}

void groupClusters(vector<size_t> &clusters)
{
    unordered_map<size_t, vector<size_t>> clusters_grp;
    size_t i = 0;
    for (size_t clus : clusters)
    {
        clusters_grp[clus].push_back(i);
        i++;
    }
    for (const auto &p : clusters_grp)
    {
        std::cerr << "*" << p.first << "\n";
        for (auto &f : p.second)
        {
            std::cerr << f << '\n';
        }
    }
}

size_t getSingleton(vector<size_t> &clusters)
{
    unordered_map<size_t, size_t> remap;
    for (size_t clus : clusters)
    {
        if (remap.count(clus))
        {
            remap[clus]++;
        }
        else
        {
            remap[clus] = 1;
        }
    }

    // for (auto &p : remap)
    //     cerr << " " << p.first << " => " << p.second << '\n';

    auto p = min_element(remap.begin(), remap.end(),
                         [](const auto &l, const auto &r)
                         { return l.second < r.second; });

    // cerr << "min " << p->second << '\n';
    // fprintf(stderr, "Output %zu clusters\n", remap.size());

    return p->second;
}

template <typename tree_type, typename signature_type>
vector<size_t> clusterSignatures(const vector<signature_type> &seqs, size_t seqCount)
{
    vector<size_t> clusters(cap);

    size_t firstNodes = 1;
    if (firstNodes > seqCount)
        firstNodes = seqCount;

    size_t offset = 0;

    tree_type tree(partree_capacity, seqs);
    tree.resizeMeans();

    if (tree_meta.readTree_)
    {
        offset = readTree(tree_meta.inputTreePath, tree);

        size_t leaves = loadSeqIDs(tree_meta.inputTreePath, tree);
        fprintf(stderr, "Loaded %zu leaves (& ambis)\n", leaves);

        // tree.printTreeJson(stderr);
    }
    offset += firstNodes;
    vector<size_t> insertionList(partree_capacity - offset); // potential nodes idx except root; root is always 0
    vector<size_t> foo(cap);

// node 0 reserved for root, node 1 reserved for leaves idx
#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 0; i < partree_capacity - offset; i++)
        {
            // insertionList.push_back(partree_capacity - i);
            insertionList[i] = partree_capacity - i - 1;
        }

#pragma omp for
        for (int i = 0; i < cap; i++)
        {
            foo[i] = i;
        }
    }

    if (random_)
    {
        // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(foo.begin(), foo.end(), default_random_engine(seed));
    }

    if (!skip_)
    {
        clusters[foo[0]] = tree.first_insert(insertionList, foo[0]);
        if (force_split_)
        {
            for (size_t i = 1; i < cap; i++)
            {
                printMsg("inserting %zu\n", foo[i]);
                size_t clus = tree.insertSplitRoot(insertionList, foo[i]);
                clusters[foo[i]] = clus;
            }
        }
        else
        {
            for (size_t i = 1; i < cap; i++)
            {
                printMsg("inserting %zu\n", foo[i]);
                size_t clus = tree.insert(insertionList, foo[i]);
                clusters[foo[i]] = clus;
            }
        }

        fprintf(stderr, "Done initial construction\nUsed %zu nodes\n", insertionList.back());
        fprintf(stderr, "Ignored %zu potential sigs\n", tree.potential_sigs.size());
        tree.updateTree();

        if (tree_meta.writeTree_)
        {
            tree.printTree(tree_meta.outputTreePath, insertionList);
            string fileName = tree_meta.outputTreePath;
            fileName = fileName + "clusters.txt";
            FILE *cFile = fopen(fileName.c_str(), "w");
            tree.printLeaves(cFile);
            fprintf(stderr, "\n\nPrinted initial tree\n");
        }
    }

    if (iteration_given)
    {
        size_t temp_singleton = singleton;
        singleton = 0;
        tree.removeAmbi(insertionList);
        // tree.removeRedundant();
        singleton = temp_singleton;
        printMsg("\n\nReinserting ambi (all)\n");
        if (updateAll_)
        {
            tree.updateAll();
        }
        else
        {
            tree.prepReinsert();
        }

#pragma omp parallel for
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = tree.reinsert(foo[i]);
            printMsg("\n Reinsert %zu at %zu\n", foo[i], clus);
            clusters[foo[i]] = clus;
        }
    }

    for (size_t run = 0; run < iteration; run++)
    {
        if (debug_)
        {
            singleton++;
        }
        fprintf(stderr, "Iteration %zu (singleton = %zu)\n", run, singleton);

        tree.removeAmbi(insertionList);
        // tree.removeRedundant();
        if (updateAll_)
        {
            tree.updateAll();
        }
        else
        {
            tree.prepReinsert();
        }

#pragma omp parallel for
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = tree.reinsert(foo[i]);
            printMsg("\n found %zu at %zu\n", foo[i], clus);
            clusters[foo[i]] = clus;
        }
    }

    if (tree_meta.writeTree_ | tree_meta.writeJSON)
    {
        tree.updateTree();
        tree.printTreeJson(stdout);
    }

    // Recursively destroy all locks
    tree.destroyLocks();
    return clusters;
}

#endif