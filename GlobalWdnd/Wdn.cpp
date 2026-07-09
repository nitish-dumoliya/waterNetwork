//
//    Minotaur -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2025 The Minotaur Team.
//

/**
 * \file Wdn.cpp
 * \brief The Wdn class for solving Water Distribution Network Design Problem
 * instances using ampl (.dat) format, with an arc-reversal + neighborhood-
 * search primal heuristic.
 * \author Nitish Kumar Dumoliya, IIT Bombay
 */

#include "Wdn.h"

#include "AMPLInterface.h"
#include "Engine.h"
#include "EngineFactory.h"
#include "Environment.h"
#include "FilterSQPEngine.h"
#include "MinotaurConfig.h"
#include "Function.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "Types.h"
#include "Variable.h"
#include "Constraint.h"
#include "Objective.h"
#include "LinearFunction.h"
#include "QuadraticFunction.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Option.h"
#include "LPEngine.h"
#include "QPEngine.h"
#include "NLPEngine.h"
#include "IpoptEngine.h"
#include "OsiLPEngine.h"
#include "CGraph.h"
#include "CNode.h"
#include "OpCode.h"
#include "Constraint.h"
#include "Solution.h"
#include "Solver.h"
#include "TreeManager.h"
#include "SOS1Handler.h"
#include "SOS2Handler.h"

#include "LinFeasPump.h"
#include "LinearHandler.h"
#include "Logger.h"

#include "Presolver.h"
#include "LinearHandler.h"
#include "NlPresHandler.h"

#include "AMPLHessian.h"
#include "AMPLInterface.h"
#include "AMPLJacobian.h"


#include <iomanip>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <algorithm>
#include <random>
#include <vector>


using namespace Minotaur;

const std::string Wdn::me_ = "Wdn: ";

namespace {
const double kFlowEps    = 1e-4;   // treat |q| below this as ~zero flow
const double kImproveTol = 1e-4;   // required objective drop to accept a move
}

Wdn::Wdn(EnvPtr env)
  : Qmax_(0.0),
    objSense_(1.0),
    status_(NotStarted),
    sol_(NULL),
    currentCost_(std::numeric_limits<double>::infinity()),
    startTime_(0.0),
    bestHitTime_(0.0),
    timeLimit_(600.0),
    numberOfNlp_(0),
    maxFailStreak_(5),
    etaLMin_(0.2),
    etaHMin_(0.2),
    alphaQMin_(0.2),
    etaLExpand_(1.2),
    etaHExpand_(1.2),
    alphaExpand_(1.2)
{
  env_ = env;
  iface_ = 0;
  ownIface_ = true;
}

Wdn::~Wdn()
{
}

void Wdn::doSetup()
{
  setInitialOptions_();
}

int Wdn::getEngine_(Engine **e)
{
  EngineFactory efac(env_);
  bool cont = false;
  EnginePtr eptr = 0;
  int err = 0;

  p_->calculateSize();
  if (p_->isLinear()) {
    eptr = efac.getLPEngine();
    if (!eptr) {
      cont = true;
    }
  }

  if (true == cont || p_->isQP()) {
    eptr = efac.getQPEngine();
    if (!eptr) {
      cont = true;
    }
  }

  if (!eptr) {
    eptr = efac.getNLPEngine();
  }

  if (!eptr) {
    env_->getLogger()->errStream()
        << "No engine available for this problem." << std::endl
        << "exiting without solving" << std::endl;
    err = 1;
  } else {
    env_->getLogger()->msgStream(LogExtraInfo)
        << me_ << "engine used = " << eptr->getName() << std::endl;
  }
  *e = eptr;
  return err;
}

NLPEnginePtr Wdn::getNLPEngine_()
{
  EngineFactory efac(env_);
  NLPEnginePtr e = efac.getNLPEngine();
  if (!e) {
    env_->getLogger()->errStream()
        << me_ << "Cannot find an NLP engine. Cannot proceed!" << std::endl;
  }
  std::cout << me_ << "Solver selected:" << &e << std::endl;
  return e;
}

void Wdn::loadData(const std::string &fileName)
{
  std::ifstream in(fileName);
  if (!in) {
    std::cerr << "Error: cannot open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (getline(in, line)) {
    if (line.find("param: nodes:") != std::string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        Node n;
        iss >> n.id >> n.elev >> n.demand >> n.pressureMin >> n.pressureMax;
        if (iss) {
          nodes_[n.id] = n;
        }
      }
    } else if (line.find("param: arcs:") != std::string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        Arc a;
        iss >> a.startNode >> a.endNode >> a.length >> a.vmax;
        if (iss) {
          arcs_[{a.startNode, a.endNode}] = a;
        }
      }
    } else if (line.find("param: pipes:") != std::string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        Pipe p;
        iss >> p.id >> p.diameter >> p.cost >> p.roughness;
        if (iss) {
          pipes_[p.id] = p;
        }
      }
    } else if (line.find("set Source:=") != std::string::npos) {
      std::istringstream iss(line);
      std::string tmp;
      iss >> tmp >> tmp;
      int s;
      while (iss >> s) {
        sources_.insert(s);
      }
    }
  }
  if (nodes_.empty()) {
    std::cerr << "Error: No nodes found in data file " << fileName
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (arcs_.empty()) {
    std::cerr << "Error: No arcs found in data file " << fileName
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (pipes_.empty()) {
    std::cerr << "Error: No pipes found in data file " << fileName
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sources_.empty()) {
    std::cerr << "Error: No sources found in data file " << fileName
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Wdn::printData() const
{
  std::cout << "Nodes:\n";
  for (NodeConstIter it = nodes_.cbegin(); it != nodes_.cend(); ++it) {
    std::cout << "  Node=" << it->first << " elev=" << it->second.elev
              << " demand=" << it->second.demand
              << " pmin=" << it->second.pressureMin
              << " pmax=" << it->second.pressureMax << "\n";
  }
  std::cout << "\nArcs:\n";
  for (ArcConstIter it = arcs_.cbegin(); it != arcs_.cend(); ++it) {
    std::cout << "  Arc = " << it->second.startNode << " -> "
              << it->second.endNode << " length=" << it->second.length
              << " vmax=" << it->second.vmax << "\n";
    // << " pipeId=" << a.pipeId << "\n";
  }

  std::cout << "\nPipes:\n";
  for (PipeConstIter it = pipes_.cbegin(); it != pipes_.cend(); ++it) {
    std::cout << "  id=" << it->second.id << " d=" << it->second.diameter
              << " C=" << it->second.cost << " R=" << it->second.roughness
              << "\n";
  }

  std::cout << "\nSource nodes:";
  for (SourceConstIter it = sources_.begin(); it != sources_.end(); ++it) {
    std::cout << " " << *it << std::endl;
  }
}


double Wdn::maxElevation()
{
  double maxElevation_ = -std::numeric_limits<double>::infinity();
  for (SourceIter sit = sources_.begin(); sit != sources_.end(); ++sit) {
    auto it = nodes_.find(*sit);
    if (it != nodes_.end()) {
      maxElevation_ = std::max(maxElevation_, it->second.elev);
    }
  }
  return maxElevation_;
}

double Wdn::calculateQmax()
{
  double Qmax = 0.0;
  for (NodeConstIter it = nodes_.cbegin(); it != nodes_.cend(); ++it) {
    if (sources_.find(it->second.id) == sources_.end())
      Qmax += it->second.demand;
  }
  std::cout << me_ << "Calculated Qmax = " << Qmax << std::endl;
  return Qmax;
}

Problem *Wdn::buildModel(ProblemPtr p)
{
  OptionDBPtr options = env_->getOptions();
  double Qmax = calculateQmax();
  double maxE = maxElevation();
  // Decision variables h[i]
  p_ = p;
  initParams_();          // <-- populate L_, D_, C_, R_, E_, d_, eps_, maxK_, scalars

  for (NodeIter it = nodes_.begin(); it != nodes_.end(); ++it) {
    Node &node = it->second;
    std::string vname = "h_" + std::to_string(node.id);
    VariablePtr h;
    if (sources_.find(node.id) != sources_.end()) {
      h = p_->newVariable(node.elev, node.elev, Continuous, vname, VarOrig);
    } else {
      h = p_->newVariable(node.elev + node.pressureMin, INFINITY, Continuous,
                          vname, VarOrig);
    }
    hvar_[node.id] = h;
  }
  std::cout << me_ << "Decision variable h[i] created for " << nodes_.size()
            << " nodes.\n";

  // Decision variables q[i,j]
  for (ArcIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    Arc &a = it->second;
    VariablePtr q =
        p_->newVariable(-Qmax_, Qmax_, Continuous,
                        "q_" + std::to_string(a.startNode) + "_" +
                            std::to_string(a.endNode),
                        VarOrig);
    qvar_[{a.startNode, a.endNode}] = q;
  }
  std::cout << me_ << "Decision variable q[i,j] created for " << arcs_.size()
            << " arcs.\n";

  // Decision variables l[i,j,k]
  for (ArcIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    Arc &a = it->second;
    for (PipeIter pit = pipes_.begin(); pit != pipes_.end(); ++pit) {
      int pid = pit->first;
      std::string vname = "l_" + std::to_string(a.startNode) + "_" +
                          std::to_string(a.endNode) + "_" +
                          std::to_string(pid);
      VariablePtr l =
          p_->newVariable(0.0, a.length, Continuous, vname, VarOrig);
      lvar_[{a.startNode, a.endNode, pid}] = l;
    }
  }
  std::cout << me_ << "Decision variable l[i,j,k] created for "
            << arcs_.size() << " arcs × " << pipes_.size() << " pipes.\n";
  // Objective Function
  addObjective();
  // Constraints
  addConstraints();

  if (options->findBool("display_problem")->getValue()) {
    p_->write(env_->getLogger()->msgStream(LogNone), 9);
  }
  return p_;
}

void Wdn::addObjective()
{
  /*
   * Objective function:
   *
   *   minimize total_cost =
   *       ∑_{(i,j) ∈ arcs} ∑_{k ∈ pipes} l[i,j,k] * C[k]
   *
   * where:
   *   - l[i,j,k] is the length of pipe diameter k installed on arc (i,j)
   *   - C[k]     is the cost per unit length of pipe diameter k
   */
  LinearFunctionPtr lf = new LinearFunction();

  for (ArcIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    Arc &a = it->second;
    for (PipeIter pIter = pipes_.begin(); pIter != pipes_.end(); ++pIter) {
      int pid = pIter->first;
      VariablePtr l = lvar_[{a.startNode, a.endNode, pid}];
      lf->addTerm(l, C_[pid]);          // coefficient = C[k]
    }
  }

  FunctionPtr f = (FunctionPtr) new Function(lf);
  ObjectivePtr obj = p_->newObjective(f, 0.0, Minimize, "total_cost");
  std::cout << me_ << "Objective function added to the problem.\n";
}

void Wdn::addConstraints()
{
  /*
   * Constraints:
   *   1. Flow balance at nodes
   *   2. Head-loss (Hazen-Williams) per arc  [smoothed / con2]
   *   3. Total length per arc
   */

  // 1. Flow balance
  for (NodeIter nodeIter = nodes_.begin(); nodeIter != nodes_.end();
       ++nodeIter) {
    Node &n = nodeIter->second;
    if (sources_.find(n.id) != sources_.end())
      continue;
    LinearFunctionPtr lf = new LinearFunction();
    for (ArcIter arcIter = arcs_.begin(); arcIter != arcs_.end(); ++arcIter) {
      Arc &a = arcIter->second;
      if (a.startNode == n.id)
        lf->addTerm(qvar_[arcIter->first], -1.0);
      if (a.endNode == n.id)
        lf->addTerm(qvar_[arcIter->first], 1.0);
    }
    FunctionPtr f = (FunctionPtr) new Function(lf);
    std::string cname = "con1_flow_balance_" + std::to_string(n.id);
    p_->newConstraint(f, n.demand, n.demand, cname);
  }
  std::cout << me_ << "Flow conservation constraints added to the problem for " << nodes_.size() << " nodes. \n";

  // 2. Original Head-loss constraints (Hazen-Williams) -- currently disabled.
  for (ArcIter arcIter = arcs_.begin(); arcIter != arcs_.end(); ++arcIter)
  {
      Arc& a = arcIter->second;
      CGraph* nlf = new CGraph();

      CNode* c_q   = nlf->newNode(qvar_[{a.startNode, a.endNode}]);
      CNode* c_absq = nlf->newNode(OpAbs, c_q, NULL);
      CNode* c_exp  = nlf->newNode(0.852);
      CNode* c_pow  = nlf->newNode(OpPowK, c_absq, c_exp);
      CNode* c_flow = nlf->newNode(OpMult, c_q, c_pow);
      std::vector<CNode *> pipeTerms;
      for (PipeIter pipeIter = pipes_.begin(); pipeIter != pipes_.end(); ++pipeIter) {
          Pipe &pid = pipeIter->second;
          CNode* c_l = nlf->newNode(lvar_[{a.startNode, a.endNode, pid.id}]);
          double coeff = 10.67 / (std::pow(pid.roughness, 1.852) * std::pow(pid.diameter, 4.87));
          CNode *c_coeff = nlf->newNode(coeff);
          CNode *c_term = nlf->newNode(OpMult, c_l, c_coeff);
          pipeTerms.push_back(c_term);
      }
      CNode *c_rSum = nlf->newNode(Minotaur::OpSumList, pipeTerms.data(), pipeTerms.size());
      CNode* c_loss = nlf->newNode(OpMult, c_flow, c_rSum);
      nlf->setOut(c_loss);
      nlf->finalize();

      LinearFunctionPtr headDiff = new LinearFunction(); 
      headDiff->addTerm(hvar_[a.startNode], 1.0);
      headDiff->addTerm(hvar_[a.endNode], -1.0);

      FunctionPtr totalFunc = (FunctionPtr) new Function(headDiff, nlf); 
      std::string cname = "hl_" + std::to_string(a.startNode) + "_" + std::to_string(a.endNode);
      // p_->newConstraint(totalFunc, 0.0, 0.0, cname);
  }

  // 2. Approximate Head-loss constraints (Approximate Hazen-Williams) -- con2.
  //    NOTE: we capture the ConstraintPtr per arc into con2_ so the heuristic
  //    can read its dual for the arc-reversal sensitivity score.
  for (ArcIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
    int i = ait->second.startNode;
    int j = ait->second.endNode;
    float eps = eps_[ait->first];
    CGraph *nlf = new CGraph();
    CNode *node1 = nlf->newNode(hvar_[i]);
    CNode *node2 = nlf->newNode(hvar_[j]);
    CNode *node3 = nlf->newNode(qvar_[{i, j}]);
    CNode *node4 = nlf->newNode(0.426);
    CNode *node5 = nlf->newNode(2);
    CNode *node6 = nlf->newNode(3);
    CNode *node7 = nlf->newNode(eps);
    CNode *node8 = nlf->newNode(eps * 0.426);
    CNode *node9 = nlf->newNode(OpPowK, node3, node5);
    CNode *node10 = nlf->newNode(OpPlus, node9, node7);
    CNode *node11 = nlf->newNode(OpPowK, node10, node4);
    CNode *node12 = nlf->newNode(OpPowK, node3, node6);
    CNode *node13 = nlf->newNode(OpMult, node12, node11);
    CNode *node14 = nlf->newNode(OpPlus, node9, node8);
    CNode *node15 = nlf->newNode(OpDiv, node13, node14);
    std::vector<CNode *> pipeTerms;
    for (PipeIter pit = pipes_.begin(); pit != pipes_.end(); ++pit) {
      int k = pit->first;
      Pipe &pipe = pit->second; 
      CNode *node16 = nlf->newNode(lvar_[{i, j, k}]);
      double coeff = 10.67 / (std::pow(pipe.roughness, 1.852) * std::pow(pipe.diameter, 4.87));
      CNode *node17 = nlf->newNode(coeff);
      CNode *node18 = nlf->newNode(OpMult, node16, node17);
      pipeTerms.push_back(node18);
    }
    CNode *c_pipeSum = nullptr;
    c_pipeSum = nlf->newNode(Minotaur::OpSumList, pipeTerms.data(), pipeTerms.size());
    CNode *c_diff = nlf->newNode(OpMinus, node1, node2);
    CNode *c_loss = nlf->newNode(OpMult, node15, c_pipeSum);
    CNode *c_final = nlf->newNode(OpMinus, c_diff, c_loss);

    nlf->setOut(c_final);
    nlf->finalize();
    FunctionPtr f = (FunctionPtr) new Function(nlf);
    std::string cname = "con2_hl_" + std::to_string(i) + "_" + std::to_string(j);
    con2_[{i, j}] = p_->newConstraint(f, 0.0, 0.0, cname);   // <-- store handle
  }
  std::cout << me_ << "Head loss constraints added to the problem for " << arcs_.size() << " arcs.\n";  

  // 3. Total length
  for (ArcIter arcIter = arcs_.begin(); arcIter != arcs_.end(); ++arcIter) {
    LinearFunctionPtr lf = new LinearFunction();
    Arc &a = arcIter->second;
    for (PipeIter pipeIter = pipes_.begin(); pipeIter != pipes_.end();
         ++pipeIter) {
      int pid = pipeIter->first;
      lf->addTerm(lvar_[{a.startNode, a.endNode, pid}], 1.0);
    }
    FunctionPtr f = (FunctionPtr) new Function(lf);
    std::string cname = "con3_totalLength_" + std::to_string(a.startNode) +
                        "_" + std::to_string(a.endNode);
    p_->newConstraint(f, a.length, a.length, cname);
  }
  std::cout << me_ << "Total length constraints added to the problem for " << arcs_.size() << " arcs. \n";
} 

double Wdn::getWallTime()
{
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}


void Wdn::initParams_()
{
  // Per-node params: E[i], D[i], P[i]
  for (NodeConstIter it = nodes_.begin(); it != nodes_.end(); ++it) {
    int i = it->second.id;
    E_[i] = it->second.elev;
    D_[i] = it->second.demand;
    P_[i] = it->second.pressureMin;   // required pressure head
  }

  // Per-arc params: L[i,j]
  for (ArcConstIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    L_[it->first] = it->second.length;
  }

  // Per-pipe params: d[k], C[k], R[k]
  for (PipeConstIter it = pipes_.begin(); it != pipes_.end(); ++it) {
    int k = it->second.id;
    d_[k] = it->second.diameter;
    C_[k] = it->second.cost;
    R_[k] = it->second.roughness;
  }

  // Scalar reductions over non-source nodes: D_min, D_max
  Dmin_ =  std::numeric_limits<double>::infinity();
  Dmax_ = -std::numeric_limits<double>::infinity();
  for (NodeConstIter it = nodes_.begin(); it != nodes_.end(); ++it) {
    if (sources_.count(it->second.id)) continue;
    Dmin_ = std::min(Dmin_, it->second.demand);
    Dmax_ = std::max(Dmax_, it->second.demand);
  }

  // Scalar reductions over pipes: d_min, d_max, R_min, R_max
  dmin_ = Rmin_ =  std::numeric_limits<double>::infinity();
  dmax_ = Rmax_ = -std::numeric_limits<double>::infinity();
  for (PipeConstIter it = pipes_.begin(); it != pipes_.end(); ++it) {
    dmin_ = std::min(dmin_, it->second.diameter);
    dmax_ = std::max(dmax_, it->second.diameter);
    Rmin_ = std::min(Rmin_, it->second.roughness);
    Rmax_ = std::max(Rmax_, it->second.roughness);
  }

  // Q_max (scalar): sum of demands over non-source nodes
  Qmax_ = calculateQmax();

  // Per-arc MaxK and eps
  const double omega = 10.67;   // set to your omega if different
  const double denom = std::pow(Rmin_, 1.852) * std::pow(dmin_, 4.87);
  for (ArcConstIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    double mk = omega * it->second.length / denom;   // MaxK[i,j]
    // maxK_[it->first] = mk;
    eps_[it->first]  = 0.0535 * std::pow(1e-3 / mk, 0.54);   // eps[i,j]
  }
}

// Adopt the engine's current solution as the incumbent if it is usable and
// strictly better than the best cost so far. Returns true if the incumbent
// was updated.
bool Wdn::updateIncumbent_(EnginePtr e)
{
  EngineStatus st = e->solve();   // or track the return of solve()
  if (st != ProvenLocalOptimal && st != ProvenOptimal) {
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "solve not optimal (" << e->getStatusString()
        << "); incumbent unchanged\n";
    return false;
  }

  double cost = e->getSolutionValue();
  if (cost >= currentCost_ - kImproveTol) {
    return false;                 // no strict improvement
  }

  // Better point found: copy primal into qval_/hval_/lval_ and record it.
  readSolution_(e);
  currentCost_ = cost;
  bestHitTime_ = getWallTime() - startTime_;

  // Keep a SolutionPtr copy so getSol()/getUb()/writeSol see the incumbent.
  ConstSolutionPtr cs = e->getSolution();
  sol_ = (SolutionPtr) new Solution(cs);

  env_->getLogger()->msgStream(LogInfo)
      << me_ << "incumbent updated: cost = " << currentCost_
      << "  time = " << bestHitTime_ << "s\n";
  return true;
}


// Write the Ipopt options file that IpoptEngine's application reads on init.
// void Wdn::writeIpoptOpt_()
// {
//   std::ofstream opt("ipopt.opt");
//   if (!opt) {
//     env_->getLogger()->errStream()
//         << me_ << "could not write ipopt.opt\n";
//     return;
//   }
//   opt << "mu_init 1e-5\n"
//       << "bound_push 0.1\n"
//       << "bound_frac 0.1\n"
//       << "linear_solver ma57\n"
//       << "ma57_pivot_order 2\n"
//       << "ma57_pivtol 1e-8\n"
//       << "ma57_pivtolmax 1e-4\n"
//       << "ma57_automatic_scaling yes\n"
//       << "ma57_pre_alloc 5\n"
//       << "mu_strategy adaptive\n"
//       << "alpha_for_y_tol 1e-6\n"
//       << "bound_relax_factor 0\n"
//       << "corrector_type affine\n"
//       << "skip_corr_if_neg_curv no\n"
//       << "hessian_approximation exact\n"
//       << "nlp_scaling_method gradient-based\n"
//       << "max_iter 2000\n"
//       << "print_user_options yes\n";   // echoes resolved options so you can verify
//   opt.close();
//   env_->getLogger()->msgStream(LogInfo)
//       << me_ << "wrote ipopt.opt\n";
// }

// Evaluate the incumbent (qval_/hval_/lval_) against the TRUE constraints
// (exact Hazen-Williams, eps = 0). Prints a per-type report and returns the
// maximum absolute violation. A small value => the incumbent is a valid UB.
double Wdn::checkConstraintViolations_()
{
  double maxFlow = 0.0, maxHl = 0.0, maxLen = 0.0, maxBnd = 0.0;

  // --- 1. Flow balance: (inflow - outflow) - demand == 0  at non-source ---
  for (NodeConstIter nit = nodes_.begin(); nit != nodes_.end(); ++nit) {
    int n = nit->second.id;
    if (sources_.count(n)) continue;
    double bal = 0.0;
    for (ArcConstIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
      int i = ait->second.startNode, j = ait->second.endNode;
      double q = qval_[ait->first];
      if (j == n) bal += q;    // inflow
      if (i == n) bal -= q;    // outflow
    }
    maxFlow = std::max(maxFlow, std::fabs(bal - D_[n]));
  }

  // --- 2. True head-loss:  (h_i - h_j) - q*|q|^0.852 * sum_k coeff_k l == 0 --
  for (ArcConstIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
    int i = ait->second.startNode, j = ait->second.endNode;
    double q = qval_[ait->first];
    double rSum = 0.0;
    for (PipeConstIter pit = pipes_.begin(); pit != pipes_.end(); ++pit) {
      int k = pit->first;
      double coeff = 10.67 / (std::pow(R_[k], 1.852) * std::pow(d_[k], 4.87));
      rSum += coeff * lval_[{i, j, k}];
    }
    // double hl = q * std::pow(std::fabs(q), 0.852);          // exact, eps = 0
    double hl = std::pow(q,3) * std::pow(std::pow(q,2) + std::pow(eps_[ait->first],2), 0.426) /(std::pow(q,2) + 0.426*std::pow(q,2));          // exact, eps = 0
    double resid = (hval_[i] - hval_[j]) - hl * rSum;
    maxHl = std::max(maxHl, std::fabs(resid));
  }

  // --- 3. Total length:  sum_k l[i,j,k] - L[i,j] == 0 ---------------------
  for (ArcConstIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
    int i = ait->second.startNode, j = ait->second.endNode;
    double sumL = 0.0;
    for (PipeConstIter pit = pipes_.begin(); pit != pipes_.end(); ++pit)
      sumL += lval_[{i, j, pit->first}];
    maxLen = std::max(maxLen, std::fabs(sumL - L_[{i, j}]));
  }

  // --- 4. Bound violations (pressure floor, non-negative lengths) ---------
  for (NodeConstIter nit = nodes_.begin(); nit != nodes_.end(); ++nit) {
    int n = nit->second.id;
    if (sources_.count(n)) continue;
    double floor = nit->second.elev + nit->second.pressureMin;   // E_i + P_i
    maxBnd = std::max(maxBnd, std::max(0.0, floor - hval_[n]));   // h_i >= floor
  }
  for (LVarConstIter lit = lvar_.begin(); lit != lvar_.end(); ++lit)
    maxBnd = std::max(maxBnd, std::max(0.0, -lval_[lit->first]));  // l >= 0

  double maxAll = std::max(std::max(maxFlow, maxHl), std::max(maxLen, maxBnd));

  std::ios_base::fmtflags f(std::cout.flags());
  std::cout << std::scientific << std::setprecision(3)
            << me_ << "constraint violations (final solution):\n"
            << "    flow balance : " << maxFlow << "\n"
            << "    head loss    : " << maxHl   << "\n"
            << "    total length : " << maxLen  << "\n"
            << "    bounds       : " << maxBnd  << "\n"
            << "    MAX          : " << maxAll  << "\n";
  std::cout.flags(f);
  return maxAll;
}

int Wdn::solve(ProblemPtr p)
{
  double walltime_start = getWallTime();
  // writeIpoptOpt_();          // <-- emit options before the engine is created

  EnginePtr engine = 0;
  PresolverPtr pres;
  VarVector *orig_v = 0;
  HandlerVector handlers;
  int err = 0;

  OptionDBPtr options = env_->getOptions();
  setInitialOptions_();
  p_ = p;
  p_->calculateSize();

  if (options->findBool("display_problem")->getValue()) {
    p_->write(env_->getLogger()->msgStream(LogNone), 12);
  }

  if (env_->getOptions()->findInt("log_level")->getValue() >= 3) {
    options->findBool("display_size")->setValue(true);
    options->findBool("display_presolved_size")->setValue(true);
  }

  if (options->findBool("display_size")->getValue()) {
    p_->writeSize(env_->getLogger()->msgStream(LogNone));
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "Starting constraint classification\n";
    p_->classifyCon();
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "Finished constraint classification\n";
  }

  if (p_->getObjective() &&
      p_->getObjective()->getObjectiveType() == Maximize) {
    objSense_ = -1.0;
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "objective sense: maximize (converted to Minimize)\n";
  } else {
    objSense_ = 1.0;
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "objective sense: minimize\n";
  }

  // get presolver.
  orig_v = new VarVector(p_->varsBegin(), p_->varsEnd());
  pres = presolve_(handlers);
  for (HandlerVector::iterator it = handlers.begin(); it != handlers.end();
       ++it) {
    delete (*it);
  }
  handlers.clear();
  if (Finished != pres->getStatus() && NotStarted != pres->getStatus()) {
    env_->getLogger()->msgStream(LogInfo)
        << me_
        << "status of presolve: " << getSolveStatusString(pres->getStatus())
        << std::endl;
    if (pres->getSolution()) {
      sol_ = pres->getPostSol(pres->getSolution());
    }
    writeSol(env_, orig_v, pres, sol_, pres->getStatus(), iface_);
    goto CLEANUP;
  }
  err = getEngine_(&engine);
  if (err) {
    goto CLEANUP;
  }

  p_->setNativeDer();
  engine->load(p_);

  startTime_   = getWallTime();
  currentCost_ = std::numeric_limits<double>::infinity();

  engine->solve();
  updateIncumbent_(engine);      // <-- adopt the first solution as incumbent

  displaySolution(engine);

  // ---- Primal heuristic: arc reversal + neighborhood search --------------
  runHeuristic_(engine);

  {
  double viol = checkConstraintViolations_();
  const double feasTol = 1e-4;
  if (viol > feasTol) {
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "WARNING: final solution violates true constraints by "
        << viol << " (> " << feasTol << "); UB may not be valid\n";
  } else {
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "final solution is feasible for the true model; "
        << "cost = " << currentCost_ << " is a valid upper bound\n";
  }
  }

  displaySolution(engine);

CLEANUP:
  if (engine)
    delete engine;
  if (pres)
    delete pres;
  if (orig_v)
    delete orig_v;

  env_->getLogger()->msgStream(LogInfo)
      << me_
      << "wall clock time used (s) = " << getWallTime() - walltime_start
      << std::endl;

  return 0;
}


void Wdn::displaySolution(EnginePtr /*e*/)
{
  std::ios_base::fmtflags f(std::cout.flags());
  std::streamsize p = std::cout.precision();

  std::cout << "Objective value = "
            << std::fixed << std::setprecision(4) << currentCost_ << std::endl;

  std::cout.flags(f);
  std::cout.precision(p);   // restore cout's original formatting
}

// void Wdn::displaySolution(EnginePtr e)
// {
//   // std::cout << "Status: " << e->getStatusString() << std::endl;
//   std::cout << "Objective value = " << currentCost_ << std::endl;
//   // std::cout << "Objective value = " << e->getSolutionValue() << std::endl;
// }

SolveStatus Wdn::getStatus()
{
  return status_;
}

PresolverPtr Wdn::presolve_(HandlerVector &handlers)
{
  PresolverPtr pres = 0;  // NULL

  p_->calculateSize();
  if (env_->getOptions()->findBool("presolve")->getValue() == true) {
    LinHandlerPtr lhandler = (LinHandlerPtr) new LinearHandler(env_, p_);
    handlers.push_back(lhandler);
    if (p_->isQP() || p_->isQuadratic() || p_->isLinear() ||
        true ==
            env_->getOptions()->findBool("use_native_cgraph")->getValue()) {
      lhandler->setPreOptPurgeVars(true);
      lhandler->setPreOptPurgeCons(true);
      lhandler->setPreOptCoeffImp(true);
    } else {
      lhandler->setPreOptPurgeVars(false);
      lhandler->setPreOptPurgeCons(false);
      lhandler->setPreOptCoeffImp(false);
    }
    if (iface_ && iface_->getNumDefs() > 0) {
      lhandler->setPreOptDualFix(false);
    } else {
      lhandler->setPreOptDualFix(true);
    }

    if (!p_->isLinear() &&
        true ==
            env_->getOptions()->findBool("use_native_cgraph")->getValue() &&
        true == env_->getOptions()->findBool("nl_presolve")->getValue()) {
      NlPresHandlerPtr nlhand =
          (NlPresHandlerPtr) new NlPresHandler(env_, p_);
      handlers.push_back(nlhand);
    }

    env_->getLogger()->msgStream(LogExtraInfo)
        << me_ << "handlers used in presolve:" << std::endl;
    for (HandlerIterator h = handlers.begin(); h != handlers.end(); ++h) {
      env_->getLogger()->msgStream(LogExtraInfo)
          << me_ << (*h)->getName() << std::endl;
    }
  }

  pres = (PresolverPtr) new Presolver(p_, env_, handlers);
  pres->standardize();
  if (env_->getOptions()->findBool("presolve")->getValue() == true) {
    pres->solve();
    for (HandlerVector::iterator h = handlers.begin(); h != handlers.end();
         ++h) {
      (*h)->writeStats(env_->getLogger()->msgStream(LogExtraInfo));
    }
  }
  return pres;
}

void Wdn::setInitialOptions_()
{
  OptionDBPtr options = env_->getOptions();
  options->findString("interface_type")->setValue("ampl");
  options->findBool("presolve")->setValue(true);
  options->findBool("nl_presolve")->setValue(true);
  options->findBool("lin_presolve")->setValue(true);
  options->findBool("msheur")->setValue(true);
  options->findString("nlp_engine")->setValue("ipopt");
}

void Wdn::writeSol(EnvPtr env, VarVector *orig_v, PresolverPtr pres,
                   SolutionPtr sol, SolveStatus status,
                   MINOTAUR_AMPL::AMPLInterface *iface)
{
  if (sol) {
    sol = pres->getPostSol(sol);
  }

  if (env->getOptions()->findFlag("ampl")->getValue() ||
      true == env->getOptions()->findBool("write_sol_file")->getValue()) {
    iface->writeSolution(sol, status);
  } else if (sol && env->getLogger()->getMaxLevel() >= LogExtraInfo &&
             env->getOptions()->findBool("display_solution")->getValue()) {
    sol->writePrimal(env->getLogger()->msgStream(LogExtraInfo), orig_v);
  }
}

void Wdn::showHelp() const
{
  std::cout
      << "Global optimization for water distribution network design "
         "problems\n"
      << "Usage:\n"
      << "To show version: wdnd -v (or --display_version yes)\n"
      << "To show all options: wdnd -= (or --display_options yes)\n"
      << "To solve an instance: wdnd --option1 [value] --option2 [value] "
         "... .dat-file\n";
}

int Wdn::showInfo()
{
  OptionDBPtr options = env_->getOptions();

  if (options->findBool("display_options")->getValue() ||
      options->findFlag("=")->getValue()) {
    options->write(std::cout);
    return 1;
  }
  if (options->findBool("display_help")->getValue() ||
      options->findFlag("?")->getValue()) {
    showHelp();
    return 1;
  }
  if (options->findBool("display_version")->getValue() ||
      options->findFlag("v")->getValue()) {
    env_->getLogger()->msgStream(LogNone)
        << me_ << "Minotaur version " << env_->getVersion() << std::endl
        << me_
       << "Global optimization for water distribution network design "
           "problems\n";
    return 1;
  }
  if (options->findString("problem_file")->getValue() == "") {
    showHelp();
    return 1;
  }
  env_->getLogger()->msgStream(LogInfo)
      << me_ << "Minotaur version " << env_->getVersion() << std::endl
      << me_
      << "Global optimization for water distribution network design "
         "problems\n";
  return 0;
}



// =============================================================================
//                         PRIMAL HEURISTIC
//   Arc reversal (iterate_acyclic_flows) + adaptive variable neighborhood search
// =============================================================================

void Wdn::readSolution_(EnginePtr e)
{
  ConstSolutionPtr s = e->getSolution();
  if (!s) return;                       // no solution to read

  const double *x = s->getPrimal();
  if (!x) return;                       // primal not available

  for (QVarConstIter it = qvar_.begin(); it != qvar_.end(); ++it)
    qval_[it->first] = x[it->second->getIndex()];
  for (HVarConstIter it = hvar_.begin(); it != hvar_.end(); ++it)
    hval_[it->first] = x[it->second->getIndex()];
  for (LVarConstIter it = lvar_.begin(); it != lvar_.end(); ++it)
    lval_[it->first] = x[it->second->getIndex()];
}

void Wdn::setInitialPoint_()
{
  std::vector<double> x0(p_->getNumVars(), 0.0);
  for (QVarConstIter it = qvar_.begin(); it != qvar_.end(); ++it)
    x0[it->second->getIndex()] = qval_[it->first];
  for (HVarConstIter it = hvar_.begin(); it != hvar_.end(); ++it)
    x0[it->second->getIndex()] = hval_[it->first];
  for (LVarConstIter it = lvar_.begin(); it != lvar_.end(); ++it)
    x0[it->second->getIndex()] = lval_[it->first];
  p_->setInitialPoint(&x0[0]);
}

double Wdn::solveNLP_(EnginePtr e)
{
  // Re-read bounds + initial point after any modification, then solve.
  e->load(p_);
  EngineStatus st = e->solve();
  ++numberOfNlp_;
  if (st != ProvenLocalOptimal && st != ProvenOptimal) {
    return std::numeric_limits<double>::infinity();
  }
  return e->getSolutionValue();
}

double Wdn::medianAbsFlow_()
{
  std::vector<double> af;
  for (ArcConstIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    double v = std::fabs(qval_[it->first]);
    if (v > 1e-6) af.push_back(v);
  }
  if (af.empty()) return 0.0;
  std::sort(af.begin(), af.end());
  return af[af.size() / 2];
}

void Wdn::arcReversal_(EnginePtr e)
{
  // Refresh the incumbent solve so the duals below match the incumbent point.
  setInitialPoint_();
  double base = solveNLP_(e);
  if (base < currentCost_ - kImproveTol) {
    readSolution_(e);
    currentCost_ = base;
  }

  // 1. In-degree of each node under the current flow orientation.
  std::map<int, int> indeg;
  for (ArcConstIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    int i = it->second.startNode, j = it->second.endNode;
    int dest = (qval_[it->first] >= 0.0) ? j : i;
    ++indeg[dest];
  }

  // 2. Candidate arcs: destination in-degree >= 2, not fixed, not visited.
  std::vector<std::pair<int, int> > cand;
  for (ArcConstIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    int i = it->second.startNode, j = it->second.endNode;
    int dest = (qval_[it->first] >= 0.0) ? j : i;
    if (indeg[dest] < 2) continue;
    if (fixArcSet_.count(it->first)) continue;
    if (visitedArcReverse_.count(it->first)) continue;
    cand.push_back(it->first);
  }

  // 3. Sensitivity score using con2 duals: sen(i,j) = -dual(i,j)*|h[i]-h[j]|.
  ConstSolutionPtr s = e->getSolution();
  const double *dcon = s->getDualOfCons();
  std::map<std::pair<int, int>, double> sen;
  for (size_t a = 0; a < cand.size(); ++a) {
    std::pair<int, int> arc = cand[a];
    int i = arc.first, j = arc.second;
    double dual = 0.0;
    std::map<std::pair<int, int>, ConstraintPtr>::iterator cit = con2_.find(arc);
    if (dcon && cit != con2_.end())
      dual = dcon[cit->second->getIndex()];
    sen[arc] = -dual * std::fabs(hval_[i] - hval_[j]);
  }

  // 4. Rank candidates by |sensitivity|, descending.
  std::sort(cand.begin(), cand.end(),
            [&](const std::pair<int, int> &a, const std::pair<int, int> &b) {
              return std::fabs(sen[a]) > std::fabs(sen[b]);
            });

  env_->getLogger()->msgStream(LogInfo)
      << me_ << "arc-reversal: " << cand.size() << " candidate arcs\n";

  // 5. Try each candidate: flip its q-bound, re-solve, accept if improved.
  for (size_t a = 0; a < cand.size(); ++a) {
    if (getWallTime() - startTime_ >= timeLimit_) {
      env_->getLogger()->msgStream(LogInfo) << me_ << "time limit reached\n";
      return;
    }
    std::pair<int, int> arc = cand[a];
    visitedArcReverse_.insert(arc);

    VariablePtr qv = qvar_[arc];
    double oldLb = qv->getLb();
    double oldUb = qv->getUb();
    double qcur  = qval_[arc];

    // Force reversal by tightening one bound to zero (restored afterwards).
    if (qcur >= 0.0) p_->changeBound(qv, oldLb, 0.0);  // force q <= 0
    else             p_->changeBound(qv, 0.0, oldUb);  // force q >= 0

    setInitialPoint_();
    double trialCost = solveNLP_(e);

    // Read the reversed solution BEFORE restoring bounds (getIndex unaffected).
    bool improved = (trialCost < currentCost_ - kImproveTol);
    if (improved) {
      readSolution_(e);
      currentCost_ = trialCost;
      bestHitTime_ = getWallTime() - startTime_;
    }

    // Restore original bounds; an accepted incumbent already has the new sign.
    p_->changeBound(qv, oldLb, oldUb);

    if (improved) {
      env_->getLogger()->msgStream(LogInfo)
          << me_ << "  reversed (" << arc.first << "," << arc.second
          << ")  new cost = " << currentCost_ << "\n";
      neighborhoodSearch_(e);   // intensify around the improved incumbent
      arcReversal_(e);          // restart from the new incumbent
      return;
    }
  }
}

void Wdn::neighborhoodSearch_(EnginePtr e)
{
  static std::mt19937 rng(12345);
  std::uniform_real_distribution<double> U(-1.0, 1.0);

  double etaL = etaLMin_, etaH = etaHMin_, alphaQ = alphaQMin_;
  int failStreak = 0;

  while (true) {
    if (getWallTime() - startTime_ >= timeLimit_) return;

    double median = medianAbsFlow_();
    double Delta = alphaQ * median;

    // Build a perturbed warm start from the incumbent WITHOUT corrupting the
    // incumbent maps (only readSolution_ on acceptance overwrites them).
    std::vector<double> x0(p_->getNumVars(), 0.0);
    for (LVarConstIter it = lvar_.begin(); it != lvar_.end(); ++it)
      x0[it->second->getIndex()] =
          std::max(0.0, lval_[it->first] * (1.0 + etaL * U(rng)));
    for (HVarConstIter it = hvar_.begin(); it != hvar_.end(); ++it)
      x0[it->second->getIndex()] = hval_[it->first] * (1.0 + etaH * U(rng));
    for (QVarConstIter it = qvar_.begin(); it != qvar_.end(); ++it)
      x0[it->second->getIndex()] = qval_[it->first] + Delta * U(rng);
    p_->setInitialPoint(&x0[0]);

    double trialCost = solveNLP_(e);
    bool improved = (trialCost < currentCost_ - kImproveTol);

    if (improved) {
      readSolution_(e);
      currentCost_ = trialCost;
      bestHitTime_ = getWallTime() - startTime_;
      env_->getLogger()->msgStream(LogInfo)
          << me_ << "  local improvement, cost = " << currentCost_ << "\n";
      etaL = etaLMin_; etaH = etaHMin_; alphaQ = alphaQMin_;  // intensify
      failStreak = 0;
    } else {
      etaL *= etaLExpand_; etaH *= etaHExpand_; alphaQ *= alphaExpand_; // diversify
      ++failStreak;

      // Terminate when perturbation can no longer push any flow past Qmax,
      // or the fail streak is exhausted.
      bool exhausted = true;
      for (ArcConstIter it = arcs_.begin(); it != arcs_.end(); ++it) {
        if (std::fabs(qval_[it->first]) + Delta <= Qmax_) {
          exhausted = false;
          break;
        }
      }
      if (exhausted || failStreak >= maxFailStreak_) {
        env_->getLogger()->msgStream(LogInfo)
            << me_ << "  local search exits\n";
        return;
      }
    }
  }
}

int Wdn::runHeuristic_(EnginePtr e)
{
  // startTime_   = getWallTime();
  numberOfNlp_ = 0;
  Qmax_        = calculateQmax();

  // Incumbent = the first NLP already solved by solve().
  readSolution_(e);
  // currentCost_ = e->getSolutionValue();
  env_->getLogger()->msgStream(LogInfo)
      << me_ << "heuristic start, initial cost = " << currentCost_ << "\n";

  neighborhoodSearch_(e);   // improve the starting point
  arcReversal_(e);          // explore flow-direction neighborhood (recurses)

  env_->getLogger()->msgStream(LogInfo)
      << me_ << "heuristic done. best cost = " << currentCost_
      << "  NLPs = " << numberOfNlp_
      << "  best-hit time = " << bestHitTime_ << "s\n";
  return 0;
}

