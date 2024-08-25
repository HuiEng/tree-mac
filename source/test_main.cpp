#include <random>
#include "tree_main_cmdline.hpp"
#include "primary_tree.hpp"
#include "cluster.hpp"
// #include "stats.hpp"

typedef primary_tree primary_tree_type;
using namespace std;
static tree_main_cmdline args; // Command line switches and arguments

// output signature size, if file too big => return empty "sigs", use readSignaturesBatch
template <typename tree_type>
size_t testTree(const string folder)
{
    // read signatureSize
    string line;
    ifstream listStream((folder + "tree.txt").c_str());
    getline(listStream, line);
    sscanf(line.c_str(), "%zu", &signatureSize);
    signatureWidth = signatureSize * sizeof(cell_type);
    getline(listStream, line);
    sscanf(line.c_str(), "%zu", &partree_capacity);
    listStream.close();

    // readTree
    tree_type tree(partree_capacity);
    tree.resizeMeans();
    size_t offset = readTree(folder, tree);
    // fprintf(stderr, "Loaded %zu nodes\n", offset);
    // size_t leaves = loadSeqIDs(folder, tree);
    // fprintf(stderr, "Loaded %zu leaves\n", leaves);
    tree.printTreeJson(stdout);
    return offset;
}

// -i inputsigs, -q order.list, -w outputsigs
template <typename cmdline_type>
void reorderedSigs(cmdline_type args)
{
    string inputFile = args.input_arg;

    vector<cell_type> seqs;
    signatureSize = readSignatures(inputFile, seqs);
    signatureWidth = signatureSize * sizeof(cell_type);
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(stderr, "Loaded %zu seqs...signatureSize %zu\n", seqCount, signatureSize);

    const string outputPath = args.topology_out_arg;
    ofstream wf(outputPath, ios::out | ios::binary);
    writeInt(wf, signatureSize);
    size_t n = 0;

    if (args.force_split_arg)
    {
        vector<size_t> foo(seqCount);
#pragma omp parallel for
        for (size_t i = 0; i < seqCount; i++)
        {
            foo[i] = i;
        }

        shuffle(foo.begin(), foo.end(), default_random_engine(seed));
        FILE *pFile = fopen(args.topology_in_arg, "w");
        for (size_t i : foo)
        {
            fprintf(pFile, "%zu\n", i);
            wf.write((char *)&seqs[i * signatureSize], signatureSize);
        }
        n = seqCount;
    }
    else
    {
        const string inputPath = args.topology_in_arg;
        string line;
        size_t i = 0;
        ifstream listStream(inputPath);
        while (getline(listStream, line))
        {
            sscanf(line.c_str(), "%zu", &i);
            wf.write((char *)&seqs[i * signatureSize], signatureSize);
            n++;
        }
    }

    wf.close();
    fprintf(stderr, "output %zu seqs\n", n);

    // // verify
    // vector<cell_type> temp;
    // signatureSize = readSignatures(outputPath, temp);
    // signatureWidth = signatureSize * sizeof(cell_type);
    // seqCount = seqs.size() / signatureSize;
    // fprintf(stderr, "Loaded %zu temp seqs...signatureSize %zu\n", seqCount, signatureSize);

    // for (size_t n = 0; n < args.sizeCap_arg; n++)
    // {
    //     fprintf(stderr, "%zu,%.2f\n", n, calcSimilarity(&seqs[n * signatureSize], &temp[n * signatureSize], signatureSize));
    // }
}

// read a tree and insert one seq at a time to check similarity status
primary_tree_type readPrimaryTree()
{
    vector<cell_type> seqs;
    signatureSize = readSignatures(args.input_arg, seqs);
    signatureWidth = signatureSize * sizeof(cell_type);
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(stderr, "Loaded %zu seqs... signatureSize %zu\n", seqCount, signatureSize);
    if (cap == 0)
    {
        cap = seqCount;
    }

    primary_tree_type tree(partree_capacity, seqs);
    tree.resizeMeans();
    readTree(tree_meta.inputTreePath, tree);
    size_t leaves = loadSeqIDs(tree_meta.inputTreePath, tree);
    fprintf(stderr, "Loaded %zu leaves (& ambis)\n", leaves);
    return tree;
}

void testTreeStatus()
{
    primary_tree_type tree = readPrimaryTree();
    vector<size_t> foo(cap);

    if (singleton != 1)
    {
        cap = 1;
        foo.resize(cap);
        foo[0] = singleton;
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < cap; i++)
        {
            foo[i] = i;
        }

        if (random_)
        {
            // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            shuffle(foo.begin(), foo.end(), default_random_engine(seed));
        }
    }

    vector<size_t> flags(cap);
#pragma omp parallel for
    for (size_t i = 0; i < cap; i++)
    {
        printMsg("inserting %zu\n", foo[i]);
        flags[foo[i]] = tree.getStatus(foo[i]);
    }

    printStatusFlag(flags);

    // vector<size_t> children = tree.childLinks[0];
    // size_t first_child = children[0];
    
    // for (size_t child:children){
    //     fprintf(stderr,"%zu,%.4f\n", child, tree.calcSimilarityWrap(tree.getMeanSig(first_child),tree.getMeanSig(child)));
    // }
}

int test_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    ios::sync_with_stdio(false); // No sync with stdio -> faster

    if (args.print_arg && !args.debug_arg)
    {
        if (args.sizeCap_given)
        {
            setArgs(args);
            primary_tree_type tree = readPrimaryTree();
            cap = args.sizeCap_arg;
            tree.printTreeJson_noRootLeaf(stdout);
            tree.updateAll();
            tree.printTreeJson_noRootLeaf(stdout);
            // tree.printTreeJson_noRootLeaf(stdout, cap);
            // vector<size_t> leaves;
            // tree.getAllLeaves(cap, leaves);
            // for (size_t leaf : leaves)
            // {
            //     fprintf(stderr, "%zu\n", leaf);
            // }
            // fprintf(stderr, "--- printed %zu leaves for %zu\n", leaves.size(), cap);
        }
        else
        {
            testTree<primary_tree_type>(args.topology_in_arg);
        }
    }
    else if (args.debug_arg)
    {
        setArgs(args);
        testTreeStatus();
    }
    else
    {

        reorderedSigs(args);
    }
    return 0;
}
