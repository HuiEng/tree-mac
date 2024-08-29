#ifndef INCLUDE_STATS_HPP
#define INCLUDE_STATS_HPP

#include <vector>
#include "bloom_filter.hpp"
#include "read.hpp"
#include "distance.hpp"

using namespace std;
size_t max_seqCount = 100;
size_t batch_size = 5;
size_t runs = 0;

class stats
{
public:
    double mean, stdev, q1, q3, max;
    bool isEmpty = true;

    stats() {}

    void update(double m, double s, double max_)
    {
        mean = m;
        stdev = s;
        q1 = mean - 0.675 * stdev;
        q3 = mean + 0.675 * stdev;
        isEmpty = false;
        max = max_;
    }

    stats(double mean, double stdev, double max_)
    {
        update(mean, stdev, max_);
    }

    stats(double mean_, double stdev_, double q1_, double q3_, double max_, size_t size)
    {
        mean = mean_ / size;
        stdev = stdev_ / size;
        q1 = q1_ / size;
        q3 = q3_ / size;
        max = max_ / size;
    }

    void printStats()
    {
        size_t col = 5;
        string s(1 + 11 * col, '-');
        fprintf(stderr, "|%-10s|%-10s|%-10s|%-10s|%-10s|\n", "mean", "stdev", "q1", "q3", "max");
        fprintf(stderr, "%s\n", s.c_str());
        fprintf(stderr, "|%-10.2f|%-10.2f|%-10.2f|%-10.2f|%-10.2f|\n", mean, stdev, q1, q3, max);
    }
};

vector<size_t> getIndices(size_t seqCount)
{
    vector<size_t> indices(seqCount);
    iota(indices.begin(), indices.end(), 0);
    if (seqCount > max_seqCount)
    {
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(indices.begin(), indices.end(), default_random_engine(seed));
        indices.resize(max_seqCount);
        // fprintf(stderr, "Random sampling %zu seqs...\n", max_seqCount);
    }

    return indices;
}

stats summarise(vector<double> const &v)
{
    stats st;
    if (v.empty())
    {
        return st;
    }

    auto const count = static_cast<double>(v.size());

    double sum = accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / count;
    double max = *max_element(v.begin(), v.end());

    vector<double> diff(count);
    // transform(v.begin(), v.end(), diff.begin(),
    //           bind2nd(minus<double>(), mean));
    double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = sqrt(sq_sum / count);
    st.update(mean, stdev, max);

    // st.printStats();
    return st;
}

stats findStatMeans(vector<stats> inputs)
{
    // fprintf(stderr, "%zu runs...\n", inputs.size());

    double mean, stdev, q1, q3, max;

    for (stats input : inputs)
    {
        mean += input.mean;
        stdev += input.stdev;
        q1 += input.q1;
        q3 += input.q3;
        max += input.max;
    }

    return stats(mean, stdev, q1, q3, max, inputs.size());
}

// template <typename funcType>
// void calcAllStatsBatch(vector<vector<seq_type>> seqs, size_t seqCount, funcType simFunc)
// {
//     vector<size_t> indices = getIndices(seqCount);
//     size_t sample_size = indices.size();
//     vector<double> data;

//     for (size_t i = 0; i < sample_size; i++)
//     {
//         size_t idx = indices[i];
//         seq_type seqA = seqs[floor(idx / batch_size)][idx % batch_size];
//         for (size_t j = i + 1; j < sample_size; j++)
//         {
//             idx = indices[j];
//             double sim = simFunc(seqA, seqs[floor(idx / batch_size)][idx % batch_size]) * 100.0;
//             data.push_back(sim);
//         }
//     }
//     summarise(data);
// }

template <typename funcType>
stats calcAllStatsBatch(vector<seq_type> seqs, funcType simFunc)
{
    size_t seqCount = seqs.size();
    vector<stats> allStats;
    for (size_t r = 0; r < runs; r++)
    {
        vector<size_t> indices = getIndices(seqCount);
        vector<double> data;

        for (size_t i = 0; i < max_seqCount; i++)
        {
            for (size_t j = i + 1; j < max_seqCount; j++)
            {
                double sim = simFunc(seqs[indices[i]], seqs[indices[j]]) * 100.0;
                data.push_back(sim);
            }
        }
        allStats.push_back(summarise(data));
    }
    stats ans = findStatMeans(allStats);
    ans.printStats();
    return ans;
}

stats calcAllStatsWindow(vector<seq_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size());
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcMatchingWindows(seqs[indices[i]], seqs[indices[j]]) * 100.0;
            data.push_back(sim);
        }
    }
    stats ans = summarise(data);
    ans.printStats();
    return ans;
}

stats calcAllStatsGlobal(vector<seq_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size());
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcJaccardGlobal(seqs[indices[i]], seqs[indices[j]]) * 100;
            // fprintf(stdout, "%zu,%zu,%.2f\n", i, j, sim);
            data.push_back(sim);
        }
    }
    stats ans = summarise(data);
    ans.printStats();
    return ans;
}

stats calcAllStatsLocal(vector<seq_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size());
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcJaccardLocal(seqs[indices[i]], seqs[indices[j]]) * 100.0;
            data.push_back(sim);
        }
    }
    stats ans = summarise(data);
    ans.printStats();
    return ans;
}

// Jaccard stats
stats calcAllStatsKmers(vector<cell_type> seqs)
{
    vector<size_t> indices = getIndices(seqs.size() / signatureSize);
    size_t sample_size = indices.size();
    vector<double> data;

    for (size_t i = 0; i < sample_size; i++)
    {
        for (size_t j = i + 1; j < sample_size; j++)
        {
            double sim = calcSimilarity(&seqs[indices[i] * signatureSize], &seqs[indices[j] * signatureSize], signatureSize) * 100;
            data.push_back(sim);
        }
    }
    stats ans = summarise(data);
    ans.printStats();
    return ans;
}

// Jaccard stats
stats calcAllStatsKmersBatch(vector<cell_type> seqs)
{
    size_t seqCount = seqs.size() / signatureSize;
    vector<stats> allStats;
    for (size_t r = 0; r < runs; r++)
    {
        vector<size_t> indices = getIndices(seqCount);
        vector<double> data;

        for (size_t i = 0; i < max_seqCount; i++)
        {
            for (size_t j = i + 1; j < max_seqCount; j++)
            {
                double sim = calcSimilarity(&seqs[indices[i] * signatureSize], &seqs[indices[j] * signatureSize], signatureSize) * 100;
                data.push_back(sim);
            }
        }
        allStats.push_back(summarise(data));
    }

    stats ans = findStatMeans(allStats);
    ans.printStats();
    return ans;
}

// // Jaccard stats
// void calcAllStatsKmersBatch(vector<cell_type> seqs, size_t seqCount)
// {
//     vector<size_t> indices = getIndices(seqCount);
//     size_t sample_size = indices.size();
//     vector<double> data;

//     for (size_t i = 0; i < sample_size; i++)
//     {
//         size_t temp = countSetBits(&seqs[indices[i] * signatureSize], signatureSize);
//         for (size_t j = i + 1; j < sample_size; j++)
//         {
//             double sim = calcSimilarity(&seqs[indices[i] * signatureSize], &seqs[indices[j] * signatureSize], signatureSize) * 100;
//             data.push_back(sim);
//         }
//     }
//     summarise(data);
// }

// void calcAllSetBits(const vector<cell_type> &sigs)
// {
//     size_t seqCount = sigs.size() / signatureSize;
//     for (size_t i = 0; i < seqCount; i++)
//     {
//         // fprintf(stderr,"%zu,%zu,%zu\n",kmerLength,windowLength,countSetBits(&sigs[i*signatureSize]));
//         fprintf(stderr, "BF density for the first seq: %zu\n", countSetBits(&sigs[i * signatureSize], signatureSize));
//         break;
//     }
// }

// stats return in percentage, coverts to decimal
double getSplitThreshold(const string bfIn, bool isList, size_t sample_size = 100, size_t runs_ = 10)
{

    vector<seq_type> seqs;
    if (isList)
    {
        seqs = readListPartitionBFSample(bfIn, signatureSize, sample_size * runs_);
    }
    else
    {
        seqs = readPartitionBFSample(bfIn, signatureSize, sample_size * runs_);
    }
    fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());
    size_t seqCount = seqs.size();
    max_seqCount = min(seqCount, sample_size);
    runs = runs_;
    // return calcAllStatsBatch(seqs, &calcJaccardGlobal).q3 / 100;
    return calcAllStatsBatch(seqs, &calcMatchingWindows).q3 / 100;
}

// no longer needed
double getSplitThresholdList(const string bfIn, size_t sample_size = 100, size_t runs_ = 10)
{

    vector<seq_type> seqs = readListPartitionBFSample(bfIn, signatureSize, sample_size * runs_);
    fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());
    size_t seqCount = seqs.size();
    max_seqCount = min(seqCount, sample_size);
    runs = runs_;
    return calcAllStatsBatch(seqs, &calcJaccardGlobal).q3 / 100;
}

// stats return in percentage, coverts to decimal
double getSplitThresholdSingle(const string bfIn, size_t sample_size = 100, size_t runs_ = 10)
{
    vector<cell_type> seqs;
    signatureSize = readSignaturesSample(bfIn, seqs, sample_size * runs_);
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(stderr, "Loaded %zu seqs...\n", seqCount);
    max_seqCount = min(seqCount, sample_size);
    runs = runs_;
    stats st = calcAllStatsKmersBatch(seqs);
    fprintf(stderr, "q3 %.2f, mean %.2f\n", st.q3, st.mean);
    return st.q3 / 100;
}

// stats return in percentage, coverts to decimal
double getSplitThresholdListSingle(const string bfIn, size_t sample_size = 100, size_t runs_ = 10)
{
    vector<cell_type> seqs;
    signatureSize = readListSample(bfIn, seqs, sample_size * runs_);
    size_t seqCount = seqs.size() / signatureSize;
    fprintf(stderr, "Loaded %zu seqs...\n", seqCount);
    max_seqCount = min(seqCount, sample_size);
    runs = runs_;
    return calcAllStatsKmersBatch(seqs).q3 / 100;
}

#endif