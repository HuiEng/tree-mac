#include <fstream>
#include <random>
#include <iostream>
#include <cstring>
// #include "minimiser.hpp"
// #include "bf_ktree.hpp"
// #include "bloom_filter.hpp"
// #include "bloom_filter3.hpp"
// #include "bf/object.hpp"

// #include "bf/all.hpp"
using namespace std;

#ifndef GIT_BRANCH
#define GIT_BRANCH "?"
#endif
#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif
string commit_branch = GIT_BRANCH;
string commit_hash = GIT_COMMIT_HASH;

static size_t kmerLength = 4;   // Kmer length
static size_t windowLength = 8; // window length

typedef int(main_func_t)(int argc, char *argv[]);
main_func_t tree_main;
main_func_t primary_tree_main;
main_func_t test_main;
main_func_t sos;

struct cmd_func
{
  const char *cmd;
  //  std::string  cmd;
  main_func_t *func;
};
cmd_func cmd_list[] = {
    {"tree", &tree_main},
    {"prim", &primary_tree_main},
    {"test", &test_main},

    /* help in all its form. Must be first non-command */
    {"help", &sos},
    {"-h", &sos},
    {"-help", &sos},
    {"--help", &sos},
    {"-?", &sos},
    {"", 0}};

void __sos(std::ostream *os)
{
  *os << "Usage: sm <cmd> [options] arg..." << std::endl
      << "Where <cmd> is one of: ";
  bool comma = false;
  for (cmd_func *ccmd = cmd_list; ccmd->func != sos; ccmd++)
  {
    *os << (comma ? ", " : "") << ccmd->cmd;
    comma = true;
  }
  *os << "." << std::endl;
  *os << "Options:" << std::endl
      << "  --help           Display this message" << std::endl;
}

int sos(int argc, char *argv[])
{
  __sos(&std::cout);
  return 0;
}

int main(int argc, char *argv[])
{

  fprintf(stderr, "current branch: %s \n", commit_branch.c_str());
  fprintf(stderr, "current commit: %s \n", commit_hash.c_str());
  std::string error;

  if (argc < 2)
  {
    error = "Too few arguments";
  }
  else
  {
    for (cmd_func *ccmd = cmd_list; ccmd->func != 0; ccmd++)
    {
      if (!strcmp(ccmd->cmd, argv[1]))
        //      if(!ccmd->cmd.compare(argv[1]))
        return ccmd->func(argc - 1, argv + 1);
    }
    error = "Unknown command '";
    error += argv[1];
    error += "'\n";
  }

  std::cerr << error << std::endl;
  __sos(&std::cerr);
  return 1;
}
