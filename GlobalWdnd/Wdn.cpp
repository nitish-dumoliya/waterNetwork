//
//    Minotaur -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2025 The Minotaur Team.
//

/**
 * \file Wdn.cpp
 * \brief The Wdn class for solving Water Distribution Network Design Problem
 * instances using ampl (.dat) format.
 * \author Ashutosh Mahajan, IIT Bombay
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


using namespace Minotaur;

const std::string Wdn::me_ = "Wdn: ";

Wdn::Wdn(EnvPtr env)
  : objSense_(1.0),
    status_(NotStarted),
    sol_(NULL)
{
  env_ = env;
  iface_ = 0;
  // p_ = p;
  // p_ = (ProblemPtr) new Problem(env_);
  // e_ = ((OsiLPEnginePtr) new OsiLPEngine(env_));
  // e_ = ((IpoptEnginePtr) new IpoptEngine(env_));
  // e_ = ((FilterSQPEnginePtr) new FilterSQPEngine(env_));
  // e_ = getNLPEngine_();
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
  // std::cout << me_ << "Calculated maximum elevation from the elevation of
  // sources = " << maxElevation_ << std::endl;
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

  for (NodeIter it = nodes_.begin(); it != nodes_.end(); ++it) {
    Node &node = it->second;
    std::string vname = "h_" + std::to_string(node.id);
    VariablePtr h;
    if (sources_.find(node.id) != sources_.end()) {
      h = p_->newVariable(node.elev, node.elev, Continuous, vname, VarOrig);
    } else {
      h = p_->newVariable(node.elev + node.pressureMin, maxE, Continuous,
                          vname, VarOrig);
    }
    hvar_[node.id] = h;
  }
  std::cout << me_ << "Decision variable h[i] created for " << nodes_.size()
            << " nodes.\n";

  // Decision variables q[i,j]
  for (ArcIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    // std::pair<int, int> key = it->first;
    Arc &a = it->second;
    VariablePtr q =
        p_->newVariable(-Qmax, Qmax, Continuous,
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
      // Pipe &p = pit->second;
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

  // if (options->findBool("display_problem")->getValue()) {
  //   p_->write(env_->getLogger()->msgStream(LogNone), 9);
  // }
  return p_;
}

void Wdn::addObjective()
{
  /*
   * Objective function:
   *
   *   minimize total_cost =
   *       ∑_{(i,j) ∈ arcs} ∑_{k ∈ pipes} l[i,j,k] * C[k] * length(i,j)
   *
   * where:
   *   - l[i,j,k] is a length of pipe diameter k is installed on arc (i,j)
   *   - C[k] is the cost per unit length of pipe diameter k
   *   - length(i,j) is the length of arc (i,j)
   *
   * This ensures that the chosen pipe design minimizes the total
   * installation cost.
   */
  ObjectivePtr obj_;
  LinearFunctionPtr lf = new LinearFunction();

  for (ArcIter it = arcs_.begin(); it != arcs_.end(); ++it) {
    for (PipeIter pIter = pipes_.begin(); pIter != pipes_.end(); ++pIter) {
      Arc &a = it->second;
      int pid = pIter->first;
      Pipe &p = pIter->second;
      VariablePtr l = lvar_[{a.startNode, a.endNode, pid}];
      lf->addTerm(l, p.cost * a.length);
    }
  }

  FunctionPtr f = (FunctionPtr) new Function(lf);
  obj_ = p_->newObjective(f, 0.0, Minimize, "total_cost");
  std::cout << me_
            << "Objective function added to the problem.\n";
}

void Wdn::addConstraints()
{
  /*
   * Constraints:
   *
   * 1. Flow balance constraint at nodes:
   *      ∑_{j: (i,j) ∈ arcs} q[i,j] - ∑_{j: (j,i) ∈ arcs} q[j,i] = demand[i]
   *
   *    Ensures conservation of flow: inflow - outflow = demand.
   *
   * 2. Head-loss (Hazen-Williams equation) constraint for each arcs:
   *      h[i] - h[j] = K * q[i,j]*|q[i,j]|^0.852
   *
   *    Where:
   *       K = ∑_{k ∈ pipes} l[i,j,k] 10.67 * l[i,j,k] / (C[k]^1.852 * d[k]^4.87)
   *
   * 3. Total length Constraint for each arcs:
   *       ∑_{k ∈ pipes} l[i,j,k] = L[i,j],   ∀ (i,j) ∈ Arcs
   *
   *     Ensure that for each arc (i,j), the sum of the selected pipe lengths
   *     across all pipe types equals the physical length of the arc.
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

  // 2. Original Head-loss constraints (Hazen-Williams)
  for (ArcIter arcIter = arcs_.begin(); arcIter != arcs_.end(); ++arcIter)
  {
      Arc& a = arcIter->second;
      CGraph* nlf = new CGraph();

      // --- Node variables (heads and flow) ---
      CNode* c_q   = nlf->newNode(qvar_[{a.startNode, a.endNode}]);  // Flow q[i,j]
      // --- Flow-dependent term: q * |q|^0.852 ---
      CNode* c_absq = nlf->newNode(OpAbs, c_q, NULL);      // |q|
      CNode* c_exp  = nlf->newNode(0.852);                 // exponent 0.852
      CNode* c_pow  = nlf->newNode(OpPowK, c_absq, c_exp); // |q|^0.852
      CNode* c_flow = nlf->newNode(OpMult, c_q, c_pow);    // q * |q|^0.852
      // --- Pipe resistance sum: Σ l[i,j,k] * coeff(k) --- //
      std::vector<CNode *> pipeTerms;
      for (PipeIter pipeIter = pipes_.begin(); pipeIter != pipes_.end(); ++pipeIter) {
          Pipe &pid = pipeIter->second;
          CNode* c_l = nlf->newNode(lvar_[{a.startNode, a.endNode, pid.id}]); // Decision variable l[i,j,k] 
          double coeff = 10.67 / (std::pow(pid.roughness, 1.852) * std::pow(pid.diameter, 4.87));
          CNode *c_coeff = nlf->newNode(coeff); // Hydraulic resistance coefficient: K = 10.67 / (R^1.852 * D^4.87)
          // Term = l[i,j,k] * coeff
          CNode *c_term = nlf->newNode(OpMult, c_l, c_coeff); // Term = l[i,j,k] * coeff
          // c_pipeSum = nlf->newNode(OpPlus, c_pipeSum, c_term);
          pipeTerms.push_back(c_term);
      }
      // --- Sum using OpSumList --- //
      CNode *c_rSum = nlf->newNode(Minotaur::OpSumList, pipeTerms.data(), pipeTerms.size());  // Σ l[i,j,k] * coeff(k)
      CNode* c_loss = nlf->newNode(OpMult, c_flow, c_rSum);
      nlf->setOut(c_loss);
      nlf->finalize();

      LinearFunctionPtr headDiff = new LinearFunction(); 

      headDiff->addTerm(hvar_[a.startNode], 1.0);  // h[i]
      headDiff->addTerm(hvar_[a.endNode], -1.0);   // h[i] - h[j]

      // h[i] - h[j] - (q * |q|^0.852)*(Σ l[i,j,k] * coeff(k))
      FunctionPtr totalFunc = (FunctionPtr) new Function(headDiff, nlf); 
      std::string cname = "hl_" + std::to_string(a.startNode) + "_" + std::to_string(a.endNode);
      // p_->newConstraint(totalFunc, 0.0, 0.0, cname);
  }

  // 2. Approximate Head-loss constraints (Approximate Hazen-Williams)
  for (ArcIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
    int i = ait->second.startNode;
    int j = ait->second.endNode;
    float eps = 1e-6;
    CGraph *nlf = new CGraph();
    // --- Node variables (heads and flow) --- //
    CNode *node1 = nlf->newNode(hvar_[i]);       // Head at node i
    CNode *node2 = nlf->newNode(hvar_[j]);       // Head at node j
    CNode *node3 = nlf->newNode(qvar_[{i, j}]);  // Flow q[i,j]
    // --- Flow-dependent term: (q^3 * (q^2 + eps)^0.426)/(q^2 + 0.426*eps) --- //
    CNode *node4 = nlf->newNode(0.426);                   // constant 0.426
    CNode *node5 = nlf->newNode(2);                       // constant 2
    CNode *node6 = nlf->newNode(3);                       // constant 3
    CNode *node7 = nlf->newNode(eps);                     // epsilon  
    CNode *node8 = nlf->newNode(eps * 0.426);             // epsilon*0.426
    CNode *node9 = nlf->newNode(OpPowK, node3, node5);    // q^2
    CNode *node10 = nlf->newNode(OpPlus, node9, node7);   // (q^2 + eps)
    CNode *node11 = nlf->newNode(OpPowK, node10, node4);  // (q^2 + eps)^0.426
    CNode *node12 = nlf->newNode(OpPowK, node3, node6);   // q^3
    CNode *node13 = nlf->newNode(OpMult, node12, node11); // q^3 * (q^2 + eps)^0.426
    CNode *node14 = nlf->newNode(OpPlus, node9, node8);   // (q^2 + 0.426*eps)
    CNode *node15 = nlf->newNode(OpDiv, node13, node14);  // (q^3 * (q^2 + eps)^0.426)/(q^2 + 0.426*eps)
    // --- Pipe resistance sum: Σ l[i,j,k] * coeff(k) --- //
    std::vector<CNode *> pipeTerms;
    for (PipeIter pit = pipes_.begin(); pit != pipes_.end(); ++pit) {
      int k = pit->first;
      Pipe &pipe = pit->second; 
      CNode *node16 = nlf->newNode(lvar_[{i, j, k}]);   // Decision variable l[i,j,k] 
      double coeff = 10.67 / (std::pow(pipe.roughness, 1.852) * std::pow(pipe.diameter, 4.87));
      CNode *node17 = nlf->newNode(coeff); // Hydraulic resistance coefficient: K = 10.67 / (R^1.852 * D^4.87)
      CNode *node18 = nlf->newNode(OpMult, node16, node17);   // Term = l[i,j,k] * coeff
      pipeTerms.push_back(node18);
    }
    // --- Sum using OpSumList --- //
    CNode *c_pipeSum = nullptr;
    c_pipeSum = nlf->newNode(Minotaur::OpSumList, pipeTerms.data(), pipeTerms.size());
    // --- Head-loss equation --- //
    CNode *c_diff = nlf->newNode(OpMinus, node1, node2);  // h_i - h_j
    CNode *c_loss = nlf->newNode(OpMult, node15, c_pipeSum);  // (q^3 * (q^2 + eps)^0.426)/(q^2 + 0.426*eps) * Σ l[i,j,k] * coeff(k)
    CNode *c_final = nlf->newNode(OpMinus, c_diff, c_loss);  // (h_i - h_j) - (q^3 * (q^2 + eps)^0.426)/(q^2 + 0.426*eps) * Σ l[i,j,k] * coeff(k)

    // --- Finalize and add constraint --- //
    nlf->setOut(c_final);
    nlf->finalize();
    FunctionPtr f = (FunctionPtr) new Function(nlf);
    std::string cname = "con2_hl_" + std::to_string(i) + "_" + std::to_string(j);
    p_->newConstraint(f, 0.0, 0.0, cname);
    // std::cout << me_ << "Constraint added for for arc(" <<i << "," << j << ")\n";
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

  // 4. Fix head at source nodes
  // for (SourceIter sIter = sources_.begin(); sIter != sources_.end();
  //      ++sIter) {
  //   Node &srcNode = nodes_.at(*sIter);
  //   LinearFunctionPtr lf = new LinearFunction();
  //   lf->addTerm(hvar_[*sIter], 1.0);
  //   FunctionPtr f = (FunctionPtr) new Function(lf);
  //   std::string cname = "con4_headFix_" + std::to_string(*sIter);
  //   p_->newConstraint(f, srcNode.elev, srcNode.elev, cname);
  // }
} 

double Wdn::getWallTime()
{
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int Wdn::solve(ProblemPtr p)
{
  double walltime_start = getWallTime();

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

  // if (false == options->findBool("use_native_cgraph")->getValue()) {
  //
  //   JacobianPtr jac = new MINOTAUR_AMPL::AMPLJacobian(iface_);
  //   p_->setJacobian(jac);
  //   HessianOfLagPtr hess = new MINOTAUR_AMPL::AMPLHessian(iface_);
  //   p_->setHessian(hess);
  //   p_->setInitialPoint(iface_->getInitialPoint(),
  //                       p_->getNumVars() - iface_->getNumDefs());
  // }

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
  engine->solve();
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

void Wdn::displaySolution(EnginePtr e)
{
  std::cout << "Status: " << e->getStatusString() << std::endl;
  std::cout << "Objective value = " << e->getSolutionValue() << std::endl;

  // for (HVarConstIter hvar_iter = hvar_.cbegin(); hvar_iter != hvar_.cend();
  // ++hvar_iter)
  //   std::cout << hvar_iter->second->getName() << " = " <<
  //   e->getSolution()->getPrimal()[hvar_iter->second->getIndex()] <<
  //   std::endl;
  //
  // for (QVarConstIter qvar_iter = qvar_.cbegin(); qvar_iter != qvar_.cend();
  // ++qvar_iter){
  //   std::cout<< qvar_iter->second->getName() << " = " <<
  //   e->getSolution()->getPrimal()[qvar_iter->second->getIndex()] <<
  //   std::endl;
  // }
  // for (LVarConstIter lvar_iter = lvar_.cbegin(); lvar_iter != lvar_.cend();
  // ++lvar_iter) {
  //   std::cout << lvar_iter->second->getName() << " = " <<
  //   e->getSolution()->getPrimal()[lvar_iter->second->getIndex()] <<
  //   std::endl;
  // }

  // return 0.0;
}

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

    // write the names.
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
