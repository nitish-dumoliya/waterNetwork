//
//    Minotaur -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2025 The Minotaur Team.
//

/**
 * \file Wdn.h
 * \brief The Wdn class for solving Water Distribution Network Design Problem instances using ampl (.dat) format.
 * \author Nitish Kumar Dumoliya, IIT Bombay
 */

#ifndef WDN_H
#define WDN_H

#include "ConBoundMod.h"
#include "Problem.h"
#include "Objective.h"
#include "NLPEngine.h"
#include "Types.h"
#include "Solver.h"
#include "Presolver.h"
#include "AMPLInterface.h"
#include "NLPEngine.h"
#include "Presolver.h"
#include "Environment.h"
#include "LPEngine.h"

#include <string>
#include <map>
#include <unordered_set>
#include <tuple>


namespace Minotaur {
/**
 * The Wdn class sets up methods for solving an Water Distribution network Design Problem instance to
 * global optimality
 */
class Wdn 
{
public:
  /// Default constructor.
  Wdn(EnvPtr env);

  /// Destroy.
  virtual ~Wdn();

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
  Problem* buildModel(ProblemPtr p);

  /// Solve optimization problem
  virtual int solve(ProblemPtr p);
  /// Calculate maximum total demand (excluding sources)
  double calculateQmax();
  /// Calculate maximum elevation from the Source nodes
  double maxElevation();
  /// get status of the last solve.
  virtual SolveStatus getStatus();

  /// get solution of the last solve.
  void displaySolution(EnginePtr e);
  virtual SolutionPtr getSol() {return sol_;};

  /// Return the upper bound for the optimal value
  double getUb();

  /// Return the lower bound for the optimal value
  double getLb();

  double getWallTime();
  
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
  double objSense_;
  // ProblemPtr oinst_;
  SolveStatus status_;
  SolutionPtr sol_;
  MINOTAUR_AMPL::AMPLInterface* iface_;
  bool ownIface_;
  // ProblemPtr inst_;
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
  
  PresolverPtr presolve_(HandlerVector &handlers); 
  int getEngine_(Engine **e);
  void writeSol(EnvPtr env, VarVector *orig_v,
            PresolverPtr pres, SolutionPtr sol, SolveStatus status,
            MINOTAUR_AMPL::AMPLInterface* iface);

};
}  // namespace Minotaur
#endif
