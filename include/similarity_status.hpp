
// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_similarity_status_HPP
#define INCLUDE_similarity_status_HPP
#include <omp.h>
#include <bitset>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "read.hpp"
#include "distance.hpp"
#include <stdarg.h>
#include <stdio.h>
#include "stats.hpp"


using namespace std;
bloom_parameters parameters;
size_t partree_capacity = 10000;
size_t singleton = 1;
size_t minClusSize = partree_capacity;
size_t tree_order = 5;
size_t calcMethod = 0;
size_t kMean_seed = 0;
bool print_ = false;
bool updateAll_ = false;

void printMsg(const char *format, ...)
{
    if (print_)
    {
        va_list args;
        va_start(args, format);
        vfprintf(stderr, format, args);
        va_end(args);
    }
}

template <typename T>
void removeVecValue(vector<T> &vec, T value)
{
    // vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());

    auto it = std::find(vec.begin(), vec.end(), value);

    // If element is found found, erase it
    if (it != vec.end())
    {
        vec.erase(it);
    }
}

template <typename T>
void removeVecIdx(vector<T> &vec, size_t idx)
{
    vec.erase(vec.begin() + idx);
}

template <typename T>
void insertVecRange(vector<T> &to, vector<T> &from)
{
    to.insert(to.end(), from.begin(), from.end());
}

template <typename T>
T returnEmpy()
{
    T empty;
    return empty;
}

enum
{
    AMBI_T = 1,
    LEAF_T = 2,
    BRANCH_T = 3,
    SUPER_T = 4,
    ROOT_T = 5
};


struct tt_data
{
    std::bitset<8> statuses;
    size_t dest = 0;
    double max_similarity = 0;
    size_t dest_branch = 0;
    double max_branch_similarity = 0;
    size_t dest_super = 0;
    double max_super_similarity = 0;
    size_t dest_root = 0;
    double max_root_similarity = 0;
    size_t mismatch = 0;

    vector<size_t> stay_branch;
    vector<size_t> stay_leaf;
    vector<size_t> nn_leaf;
    vector<size_t> nn_branch;
    vector<size_t> stay_super;
    vector<size_t> stay_root;
};



size_t dt_append(tt_data &dt_to, tt_data &dt_from)
{
    if (dt_from.statuses.none())
    {
        return 0;
    }
    dt_to.statuses |= dt_from.statuses;
    if (dt_from.max_similarity > dt_to.max_similarity)
    {
        dt_to.dest = dt_from.dest;
        dt_to.max_similarity = dt_from.max_similarity;
    }
    insertVecRange(dt_to.stay_leaf, dt_from.stay_leaf);
    insertVecRange(dt_to.nn_leaf, dt_from.nn_leaf);

    if (dt_from.max_branch_similarity > dt_to.max_branch_similarity)
    {
        dt_to.dest_branch = dt_from.dest_branch;
        dt_to.max_branch_similarity = dt_from.max_branch_similarity;
    }
    insertVecRange(dt_to.stay_branch, dt_from.stay_branch);
    insertVecRange(dt_to.nn_branch, dt_from.nn_branch);

    if (dt_from.max_super_similarity > dt_to.max_super_similarity)
    {
        dt_to.dest_super = dt_from.dest_super;
        dt_to.max_super_similarity = dt_from.max_super_similarity;
    }
    insertVecRange(dt_to.stay_super, dt_from.stay_super);
    return 1;
}

size_t dt_candidate_cnt(tt_data &dt)
{
    if (dt.statuses.none())
    {
        return 0;
    }
    return dt.stay_leaf.size() +
           dt.nn_leaf.size() +
           dt.stay_branch.size() +
           dt.nn_branch.size() +
           dt.stay_super.size();
}



size_t doBit1(tt_data dt)
{
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
        // fprintf(stdout, "2,>>Stat>> Stay Leaf-\n");
        return 2;
        break;
    case 2:
        // fprintf(stdout, "3,>>Stat>> NN Leaf-\n");
        return 3;
        break;
    case 3:
        // fprintf(stdout, "4,>>Stat>> Stay Branch-\n");
        return 4;
        break;
    case 4:
        // fprintf(stdout, "5,>>Stat>> NN Branch-%zu; mismatch:%zu\n", dt.nn_branch.size(), dt.mismatch);
        return 5;
        break;
    case 5:
        // fprintf(stdout, "6,>>Stat>> Stay Super-%zu\n", dt.stay_super.size());
        return 6;
        break;
    case 7:
        // fprintf(stdout, "7,>>Stat>> Stay Root-\n");
        return 7;
        break;
    default:
        // fprintf(stdout, "ERROR 1 bit, pos!!\n");
        return 0;
        break;
    }
}

size_t doBit2(tt_data dt)
{
    size_t tempbits = (dt.statuses & bitset<8>("00101010")).count();
    // only stay
    if (tempbits == 2)
    {
        if (dt.statuses.test(5))
        {
            if (dt.statuses.test(1))
            {
                // fprintf(stdout, "8,>>Stat>> Stay Leaf and Stay Super-\n");
                return 8;
            }
            else
            {
                // fprintf(stdout, "9,>>Stat>> Stay Branch and Stay Super-\n");
                return 9;
            }
        }
        else
        {
            // fprintf(stdout, "10,>>Stat>> Stay Branch and Stay Super-\n");
            return 10;
        }
    }
    else if (tempbits == 1)
    {
        if (dt.statuses.test(4))
        {
            if (dt.statuses.test(1))
            {
                // fprintf(stdout, "11,>>Stat>> Stay Leaf and NN Branch-\n");
                return 11;
            }
            else if (dt.statuses.test(3))
            {
                // fprintf(stdout, "12,>>Stat>> Stay Branch and NN Branch-\n");
                return 12;
            }
            else if (dt.statuses.test(5))
            { //?
              // fprintf(stdout, "13,>>Stat>> NN Branch and Stay Super-\n");
                return 13;
            }
        }
        else
        {
            if (dt.statuses.test(1))
            {
                // fprintf(stdout, "14,>>Stat>> Stay Leaf and NN Leaf-\n");
                return 14;
            }
            else if (dt.statuses.test(3))
            {
                // fprintf(stdout, "15,>>Stat>> NN Leaf and Stay Branch-\n");
                return 15;
            }
            else if (dt.statuses.test(5))
            {
                // fprintf(stdout, "16,>>Stat>> NN Leaf and Stay Super-\n");
                return 16;
            }
        }
    }
    else
    {
        // fprintf(stdout, "17,>>Stat>> NN Leaf and NN Branch-\n");
        return 17;
    }
    return 0;
}

size_t doBit3(tt_data dt)
{
    if (dt.statuses.test(5))
    {
        if (dt.statuses.test(1))
        {
            if (dt.statuses.test(2))
            {
                // fprintf(stdout, "18,>>Stat>> Stay Leaf and NN Leaf and Stay Super-\n");
                return 18;
            }
            else if (dt.statuses.test(3))
            {
                // fprintf(stdout, "19,>>Stat>> Stay Leaf and Stay Branch and Stay Super-\n");
                return 19;
            }
            else
            {
                // fprintf(stdout, "20,>>Stat>> Stay Leaf and NN branch and Stay Super-\n");
                return 20;
            }
        }
        else if (dt.statuses.test(2))
        {
            if (dt.statuses.test(3))
            {
                // fprintf(stdout, "21,>>Stat>> NN Leaf and Stay Branch and Stay Super-\n");
                return 21;
            }
            else
            {
                // fprintf(stdout, "22,>>Stat>> NN Leaf and NN branch and Stay Super-\n");
                return 23;
            }
        }
        else
        {
            // fprintf(stdout, "23,>>Stat>> Stay Branch and NN branch and Stay Super-\n");
            return 23;
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
                    // fprintf(stdout, "24,>>Stat>> Stay Leaf and NN Leaf and Stay Branch-\n");
                    return 24;
                }
                else
                {
                    // fprintf(stdout, "25,>>Stat>> Stay Leaf and NN Leaf and NN Branch-\n");
                    return 25;
                }
            }
            else
            {
                // fprintf(stdout, "26,>>Stat>> Stay Leaf and Stay Branch and NN branch-\n");
                return 26;
            }
        }
        else
        {
            // fprintf(stdout, "27,>>Stat>> NN Leaf and Stay Branch and NN branch-\n");
            return 27;
        }
    }
    return 0;
}

size_t doBit4(tt_data dt)
{
    if (dt.statuses.test(1))
    {
        if (dt.statuses.test(2))
        {
            if (dt.statuses.test(3))
            {
                if (dt.statuses.test(4))
                {
                    // fprintf(stdout, "28,>>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch-\n");
                    return 28;
                }
                else
                {
                    // fprintf(stdout, "29,>>Stat>> Stay Leaf and NN Leaf and Stay Branch and Stay Super-\n");
                    return 29;
                }
            }
            else
            {
                // fprintf(stdout, "30,>>Stat>> Stay Leaf and NN Leaf and NN Branch and Stay Super-\n");
                return 30;
            }
        }
        else
        {
            // fprintf(stdout, "31,>>Stat>> Stay Leaf and Stay Branch and NN Branch and Stay Super-\n");
            return 31;
        }
    }
    else
    {
        // fprintf(stdout, "32,>>Stat>> NN Leaf and Stay Branch and NN Branch and Stay Super-\n");
        return 32;
    }
}

size_t getStatusFlag(tt_data dt)
{
    if (dt.statuses.none())
    {
        // fprintf(stdout, "0,>>Stat>> Mismatch-\n");
        return 1;
    }

    size_t setbits = dt.statuses.count();
    if (setbits == 1)
    {
        return doBit1(dt);
    }
    else if (setbits == 2)
    {
        return doBit2(dt);
    }
    else if (setbits == 3)
    {
        return doBit3(dt);
    }
    else if (setbits == 4)
    {
        return doBit4(dt);
    }
    else if (setbits == 5)
    {
        // fprintf(stdout, "33,>>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch and Stay Super-\n");
        return 33;

        // t_branch = createBranch(node, insertionList, dt.stay_leaf);
        // for (size_t l : dt.nn_leaf)
        // {
        //     moveParent(l, t_branch);
        // }
        // fprintf(stdout, "combined leaves to branch %zu\n", t_branch);
        // dt.stay_branch.push_back(t_branch);
        // dt.dest_branch = t_branch;
        // dt.statuses.reset(1);
        // dt.statuses.reset(2);
        // doBit3(dt, signature, insertionList, idx, node);
    }
    else
    {
        // fprintf(stdout, ">>Stat>> ERROR!!\n");
        return 0;
    }
}

void printStatusFlag(vector<size_t> &flags)
{
    vector<const char *> statusInfo{">>Stat>> ERROR!!",
                                    ">>Stat>> Mismatch-",
                                    ">>Stat>> Stay Leaf-",
                                    ">>Stat>> NN Leaf-",
                                    ">>Stat>> Stay Branch-",
                                    ">>Stat>> NN Branch-",
                                    ">>Stat>> Stay Super-",
                                    ">>Stat>> Stay Root-",
                                    ">>Stat>> Stay Leaf and Stay Super-",
                                    ">>Stat>> Stay Branch and Stay Super-",
                                    ">>Stat>> Stay Branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and NN Branch-",
                                    ">>Stat>> Stay Branch and NN Branch-",
                                    ">>Stat>> NN Branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and NN Leaf-",
                                    ">>Stat>> NN Leaf and Stay Branch-",
                                    ">>Stat>> NN Leaf and Stay Super-",
                                    ">>Stat>> NN Leaf and NN Branch-",
                                    ">>Stat>> Stay Leaf and NN Leaf and Stay Super-",
                                    ">>Stat>> Stay Leaf and Stay Branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and NN branch and Stay Super-",
                                    ">>Stat>> NN Leaf and Stay Branch and Stay Super-",
                                    ">>Stat>> NN Leaf and NN branch and Stay Super-",
                                    ">>Stat>> Stay Branch and NN branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and NN Leaf and Stay Branch-",
                                    ">>Stat>> Stay Leaf and NN Leaf and NN Branch-",
                                    ">>Stat>> Stay Leaf and Stay Branch and NN branch-",
                                    ">>Stat>> NN Leaf and Stay Branch and NN branch-",
                                    ">>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch-",
                                    ">>Stat>> Stay Leaf and NN Leaf and Stay Branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and NN Leaf and NN Branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and Stay Branch and NN Branch and Stay Super-",
                                    ">>Stat>> NN Leaf and Stay Branch and NN Branch and Stay Super-",
                                    ">>Stat>> Stay Leaf and NN Leaf and Stay Branch and NN Branch and Stay Super-"};

    size_t i = 0;
    for (size_t flag : flags)
    {
        fprintf(stdout, "%zu,%zu,%s\n", i++, flag, statusInfo[flag]);
    }

    unordered_map<size_t, vector<size_t>> clusters_grp;
    i = 0;
    for (size_t clus : flags)
    {
        clusters_grp[clus].push_back(i);
        i++;
    }
    fprintf(stderr, "count,flag,info\n");
    for (const auto &p : clusters_grp)
    {
        size_t flag = p.first;
        fprintf(stderr, "%zu,%zu,%s\n", p.second.size(), flag, statusInfo[flag]);
    }
}

#endif