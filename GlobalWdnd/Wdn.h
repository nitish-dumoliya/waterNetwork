//
//    Minotaur -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2025 The Minotaur Team.
//

/**
 * \file Wdn.h
 * \brief The Wdn class for solving Water Distribution Network Design Problem
 *        instances using ampl (.dat) format, with an Arc-Reversal +
 *        Neighborhood-Search primal heuristic.
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
#include <set>
#include <vector>
#include <unordered_set>
#include <tuple>


namespace Minotaur {
/**
 * The Wdn class sets up methods for solving a Water Distribution Network Design
 * Problem instance to global optimality, and includes a primal heuristic
 * (arc reversal + adaptive neighborhood search) to produce good upper bounds.
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


  // ---- Model parameters (mirror the AMPL params; filled by initParams_) ----
  std::map<std::pair<int,int>, double> L_;      // L[i,j]   arc length
  std::map<int, double>                D_;      // D[i]     nodal demand
  std::map<int, double>                C_;      // C[k]     pipe cost
  std::map<int, double>                P_;      // P[i]     required pressure head
  std::map<int, double>                R_;      // R[k]     pipe roughness
  std::map<int, double>                E_;      // E[i]     node elevation
  std::map<int, double>                d_;      // d[k]     pipe diameter
  std::map<std::pair<int,int>, double> alpha_;  // alpha[i,j]
  std::map<std::pair<int,int>, double> eps_;    // eps[i,j]
  std::map<std::pair<int,int>, double> maxK_;   // MaxK[i,j]
  
  double Dmin_, Dmax_;   // over non-source nodes
  double dmin_, dmax_;   // over pipes
  double Rmin_, Rmax_;   // over pipes
  // (Q_max already exists as Qmax_)
  
  void initParams_();    // populate all of the above
  bool updateIncumbent_(EnginePtr e);   // adopt e's solution if it improves

  // Decision variables
  std::map<int, VariablePtr> hvar_;                  // head at node i
  std::map<std::pair<int, int>, VariablePtr> qvar_;  // flow on arc (i,j)
  std::map<std::tuple<int, int, int>, VariablePtr>
      lvar_;  // pipe selection l[i,j,k]

  // ===== Heuristic state ==================================================
  // Head-loss (con2) constraint handle per arc -- needed for dual lookup.
  std::map<std::pair<int, int>, ConstraintPtr> con2_;

  // Incumbent solution values (kept in sync with the best point found).
  std::map<std::pair<int, int>, double> qval_;        // q[i,j]
  std::map<int, double> hval_;                        // h[i]
  std::map<std::tuple<int, int, int>, double> lval_;  // l[i,j,k]

  double currentCost_;   // best (lowest) objective found so far
  double startTime_;     // wall-clock start of the heuristic
  double bestHitTime_;   // time at which the best solution was found
  double timeLimit_;     // heuristic time budget (seconds)
  int    numberOfNlp_;   // count of NLP solves
  int    maxFailStreak_; // consecutive non-improving trials before giving up

  // Adaptive VNS radii (defaults mirror the Python __init__).
  double etaLMin_;    // min pipe-length perturbation radius
  double etaHMin_;    // min head perturbation radius
  double alphaQMin_;  // min flow perturbation factor
  double etaLExpand_; // growth factor for etaL on failure
  double etaHExpand_; // growth factor for etaH on failure
  double alphaExpand_;// growth factor for alphaQ on failure

  // Arcs excluded from reversal / already tried.
  std::set<std::pair<int, int> > fixArcSet_;
  std::set<std::pair<int, int> > visitedArcReverse_;
  // ========================================================================

  // Internal methods for building the model
  void addObjective();
  void addConstraints();
  void setInitialOptions_();
  // void writeIpoptOpt_();   // emit ipopt.opt with the solver options
  double checkConstraintViolations_();   // max violation of the true constraints

  PresolverPtr presolve_(HandlerVector &handlers); 
  int getEngine_(Engine **e);
  void writeSol(EnvPtr env, VarVector *orig_v,
            PresolverPtr pres, SolutionPtr sol, SolveStatus status,
            MINOTAUR_AMPL::AMPLInterface* iface);

  // ===== Heuristic methods ================================================
  /// Copy the engine's primal solution into qval_/hval_/lval_.
  void readSolution_(EnginePtr e);
  /// Push the incumbent maps in as the NLP initial (warm-start) point.
  void setInitialPoint_();
  /// Solve p_ at current bounds/init point; return objective (or +inf).
  double solveNLP_(EnginePtr e);
  /// Median of |q| over arcs with non-trivial flow.
  double medianAbsFlow_();
  /// Arc-reversal neighborhood (port of iterate_acyclic_flows).
  void arcReversal_(EnginePtr e);
  /// Adaptive neighborhood search (port of local_solution_improvement_*).
  void neighborhoodSearch_(EnginePtr e);
  /// Top-level heuristic driver; call after the first NLP solve.
  int runHeuristic_(EnginePtr e);
  // ========================================================================

};
}  // namespace Minotaur
#endif
