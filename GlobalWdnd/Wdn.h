#ifndef WDN_H
#define WDN_H

#include "Problem.h"
#include "Objective.h"
#include "NLPEngine.h"
#include "Types.h"
#include <string>
#include <map>
#include <unordered_set>
#include <tuple>

namespace Minotaur {

  class Wdn {
  public:
    Wdn(EnvPtr env);
    ~Wdn();

    /// Starting setup for wdnd
    void doSetup();

    /// Show help messages
    void showHelp() const;

    /// Display information
    int showInfo();

    /// Load network data from file
    void loadData(const std::string &fileName);

    /// Print network data
    void printData() const;

    /// Build optimization model
    void buildModel();

    /// Solve optimization problem
    void solve();

    /// Calculate maximum total demand (excluding sources)
    double calculateQmax();

    struct Node
    {
      int id;
      double elev;
      double demand;
      double pressureMin;
      double pressureMax;
    };

    struct Arc
    {
      int startNode;  // node where arc starts
      int endNode;    // node where arc ends
      double length;
      double vmax;
      int pipeId;  // associated pipe type
    };

    struct Pipe
    {
      int id;
      double diameter;
      double cost;
      double roughness;
    };
    // ===== Iterator typedefs / aliases =====
    using NodeIter = std::map<int, Node>::iterator;
    using NodeConstIter = std::map<int, Node>::const_iterator;

    using ArcIter = std::map<std::pair<int, int>, Arc>::iterator;
    using ArcConstIter = std::map<std::pair<int, int>, Arc>::const_iterator;

    using PipeIter = std::map<int, Pipe>::iterator;
    using PipeConstIter = std::map<int, Pipe>::const_iterator;

    using SourceIter = std::unordered_set<int>::iterator;
    using SourceConstIter = std::unordered_set<int>::const_iterator;

    using HVarIter = std::map<int, VariablePtr>::iterator;
    using HVarConstIter = std::map<int, VariablePtr>::const_iterator;

    using QVarIter = std::map<std::pair<int, int>, VariablePtr>::iterator;
    using QVarConstIter =
        std::map<std::pair<int, int>, VariablePtr>::const_iterator;

    using LVarIter =
        std::map<std::tuple<int, int, int>, VariablePtr>::iterator;
    using LVarConstIter =
        std::map<std::tuple<int, int, int>, VariablePtr>::const_iterator;
    // ========================================
  private:
    const static std::string me_;
    EnvPtr env_;
    ProblemPtr p_;
    EnginePtr e_;
    NLPEnginePtr getNLPEngine_();
    double Qmax_;  // maximum total demand (excluding source)

    // Network data structures
    std::map<int, Node> nodes_;                // node ID -> Node
    std::map<std::pair<int, int>, Arc> arcs_;  // (start,end) -> Arc
    std::map<int, Pipe> pipes_;                // pipe ID -> Pipe
    std::unordered_set<int> sources_;          // node IDs that are sources

    // Decision variables
    std::map<int, VariablePtr> hvar_;                  // head at node i
    std::map<std::pair<int, int>, VariablePtr> qvar_;  // flow on arc (i,j)
    std::map<std::tuple<int, int, int>, VariablePtr>
        lvar_;  // pipe selection l[i,j,k]

    // Internal methods for building the model
    void addObjective();
    void addConstraints();
    void setInitialOptions_();
  };

}  // namespace Minotaur

#endif
