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

#include <fstream>
#include <ostream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>

using namespace Minotaur;

const std::string Wdn::me_ = "Wdnd-Global: ";

Wdn::Wdn(EnvPtr env)
  : env_(env)
{
  p_ = (ProblemPtr) new Problem(env_);
  e_ = ((OsiLPEnginePtr) new OsiLPEngine(env_));
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
    return;
  }

  std::string line;
  while (getline(in, line)) {
    if (line.find("param: nodes:") != std::string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        Node n;
        iss >> n.id >> n.elev >> n.demand >> n.pressureMin >> n.pressureMax;
        if (iss)
          nodes_[n.id] = n;  // insert into map
      }
    } else if (line.find("param: arcs:") != std::string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        Arc a;
        iss >> a.startNode >> a.endNode >> a.length >> a.vmax;
        if (iss)
          arcs_[{a.startNode, a.endNode}] = a;  // insert into map
      }
    } else if (line.find("param: pipes:") != std::string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        std::istringstream iss(line);
        Pipe p;
        iss >> p.id >> p.diameter >> p.cost >> p.roughness;
        if (iss)
          pipes_[p.id] = p;  // insert into map
      }
    } else if (line.find("set Source:=") != std::string::npos) {
      // Example: set Source := 1 4 7 ;
      std::istringstream iss(line);
      std::string tmp;
      iss >> tmp >> tmp;  // skip "set Source:="
      int s;
      while (iss >> s)
        sources_.insert(s);  // insert into unordered_set
    }
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
    std::cout << "  Arc = " << it->second.startNode << " -> " << it->second.endNode
              << " length=" << it->second.length << " vmax=" << it->second.vmax << "\n";
    // << " pipeId=" << a.pipeId << "\n";
  }

  std::cout << "\nPipes:\n";
  for (PipeConstIter it = pipes_.cbegin(); it != pipes_.cend(); ++it) {
    std::cout << "  id=" << it->second.id << " d=" << it->second.diameter << " C=" << it->second.cost
              << " R=" << it->second.roughness << "\n";
  }

  std::cout << "\nSource nodes:";
  for (SourceConstIter it = sources_.begin(); it != sources_.end(); ++it){
    std::cout << " " << *it << std::endl;
  }
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

void Wdn::buildModel()
{
  OptionDBPtr options = env_->getOptions();
  double Qmax = calculateQmax();

  // Decision variables h[i]
    for (NodeConstIter it = nodes_.cbegin(); it != nodes_.cend(); ++it) {
    std::string vname = "h_" + std::to_string(it->second.id);
    VariablePtr h = p_->newVariable(
        it->second.elev + it->second.pressureMin, it->second.elev + it->second.pressureMax, Continuous, vname);
    hvar_[it->second.id] = h;
  }
  std::cout << me_ << "Decision variable h[i] created for " << nodes_.size()
            << " nodes.\n";

  // Decision variables q[i,j]
  for (ArcConstIter it = arcs_.cbegin(); it != arcs_.cend(); ++it) {
    std::pair<int, int> key = it->first;
    // Arc &a = it->second;
    VariablePtr q = p_->newVariable(
        -Qmax, Qmax, Continuous,
        "q_" + std::to_string(it->second.startNode) + "_" + std::to_string(it->second.endNode));
    qvar_[key] = q;
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
      VariablePtr l = p_->newVariable(0.0, a.length, Continuous, vname);
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
  obj_ = p_->newObjective(f, 0.0, Minimize);
  std::cout << me_ << "Objective function successfully added to the problem.\n";
}

void Wdn::addConstraints()
{
  /*
   * Constraints:
   *
   * 1. Flow balance constraint at nodes:
   *      ∑_{j: (i,j) ∈ arcs} q[i,j] - ∑_{j: (j,i) ∈ arcs} q[j,i] =
   * demand[i]
   *
   *    Ensures conservation of flow: inflow - outflow = demand.
   *
   * 2. Head-loss (Hazen-Williams equation) constraint for each arcs:
   *      h[i] - h[j] = K * q[i,j]*|q[i,j]|^0.852
   *
   *    Where:
   *       K = ∑_{k ∈ pipes} l[i,j,k] 10.67 * l[i,j,k] / (C[k]^1.852 *
   * d[k]^4.87)
   *
   * 3. Total length Constraint for each arcs:
   *       ∑_{k ∈ pipes} l[i,j,k] = L[i,j],   ∀ (i,j) ∈ Arcs
   *
   *     Ensure that for each arc (i,j), the sum of the selected pipe lengths
   *     across all pipe types equals the physical length of the arc.
   */

  // 1. Flow balance
  for (NodeIter nodeIter = nodes_.begin(); nodeIter != nodes_.end(); ++nodeIter) {
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
    std::string cname = "flow_balance_" + std::to_string(n.id);
    p_->newConstraint(f, n.demand, n.demand, cname);
  }

  // 2. Head-loss constraints (Hazen-Williams)

  for (ArcIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
    int i = ait->second.startNode;
    int j = ait->second.endNode;
    CGraph *nlf = new CGraph();
    CNode *c_h_i, *c_h_j, *c_q, *c_absq, *c_exp, *c_pipeSum, *c_tmp; 
    c_h_i = nlf->newNode(hvar_[i]); // Node heads
    c_h_j = nlf->newNode(hvar_[j]); // Node heads 
    c_q = nlf->newNode(qvar_[{i, j}]); // Node Flow q[i,j] 
    // c_absq = nlf->newNode(OpAbs, c_q);  // Node Absolute value: |q|
    CNode *cq_sq = nlf->newNode(OpMult, c_q, c_q);  // q^2
    CNode *cq_abs = nlf->newNode(OpPowK, cq_sq, nlf->newNode(0.5));  // sqrt(q^2)
    // |q|^0.852
    c_exp = nlf->newNode(0.852);
    c_absq = nlf->newNode(OpPowK, cq_abs, c_exp);
    // q * |q|^0.852
    c_q = nlf->newNode(OpMult, c_q, c_absq);

    // Build sum over pipes
    c_pipeSum = nlf->newNode(0.0);  // initialize

    for (PipeIter pit = pipes_.begin(); pit != pipes_.end(); ++pit) {
      int k = pit->first;
      const Pipe &pipe = pit->second;
      CNode *c_l = nlf->newNode(lvar_[{i, j, k}]);  // Node for decision variable l[i,j,k] 
      double coeff = 10.67 / (std::pow(pipe.roughness, 1.852) * std::pow(pipe.diameter, 4.87)); 
      CNode *c_coeff = nlf->newNode(coeff); // Node for coefficient K = 10.67 / (R^1.852 * D^4.87) 
      c_tmp = nlf->newNode(OpMult, c_l, c_coeff); // term = l[i,j,k] * coeff
      // add to sum
      c_pipeSum = nlf->newNode(OpPlus, c_pipeSum, c_tmp);
    }
    // Multiply flow term by pipe sum
    CNode *c_loss = nlf->newNode(OpMult, c_q, c_pipeSum);
    // Constraint: (h_i - h_j) - loss = 0
    CNode *c_diff = nlf->newNode(OpMinus, c_h_i, c_h_j);
    CNode *c_final = nlf->newNode(OpMinus, c_diff, c_loss);
    nlf->setOut(c_final);
    nlf->finalize();
    FunctionPtr f = (FunctionPtr) new Function(nlf);
    std::string cname = "hl_" + std::to_string(i) + "_" + std::to_string(j);
    // p_->newConstraint(f, 0.0, 0.0, cname);
  }

  // 3. Total length
  for (ArcIter arcIter = arcs_.begin(); arcIter!=arcs_.end(); ++arcIter) {
    LinearFunctionPtr lf = new LinearFunction();
    Arc &a = arcIter->second; 
    for (PipeIter pipeIter = pipes_.begin(); pipeIter!=pipes_.end(); ++pipeIter) {
      int pid = pipeIter->first; 
      lf->addTerm(lvar_[{a.startNode, a.endNode, pid}], 1.0);
    }
    FunctionPtr f = (FunctionPtr) new Function(lf);
    std::string cname = "con3_" + std::to_string(a.startNode) + "_" + std::to_string(a.endNode);
    p_->newConstraint(f, a.length, a.length, cname);
  }

  // 4. Fix head at source nodes
  for (SourceIter sIter = sources_.begin(); sIter != sources_.end(); ++sIter) {
    const Node &srcNode = nodes_.at(*sIter);
    LinearFunctionPtr lf = new LinearFunction();
    lf->addTerm(hvar_[*sIter], 1.0);
    FunctionPtr f = (FunctionPtr) new Function(lf);
    std::string cname = "con4_headFix_" + std::to_string(*sIter);
    p_->newConstraint(f, srcNode.elev + srcNode.pressureMin, srcNode.elev + srcNode.pressureMin, cname);
  }

  // 2. Head-loss constraints(Hazen-Williams formula)
  // for (ArcIter ait = arcs_.begin(); ait != arcs_.end(); ++ait) {
  //   Arc &a = ait->second;
  //
  //   VariablePtr hi = hvar_[a.startNode];
  //   VariablePtr hj = hvar_[a.endNode];
  //   VariablePtr qij = qvar_[{a.startNode, a.endNode}];
  //   LinearFunctionPtr lf = new LinearFunction();
  //
  //   // Create nonlinear function: h_i - h_j - K * q * |q|^0.852 = 0
  //   NonlinearFunctionPtr nlf = new NonlinearFunction;
  //
  //   // Add linear terms: h_i - h_j
  //   lf->addTerm(hi, 1.0);
  //   lf->addTerm(hj, -1.0);
  //   FunctionPtr f = (FunctionPtr) new Function(lf);
  //   nlf->addTerm(hi, 1.0);
  //
  //   for (PipeIter pit = pipes_.begin(); pit != pipes_.end(); ++pit) {
  //     Pipe &p = pit->second;
  //
  //     // Get length variable
  //     VariablePtr lijk = lvar_[{a.startNode, a.endNode, p.id}];
  //
  //     // Compute Hazen-Williams coefficient
  //     double K = 10.67 * a.length / (std::pow(p.roughness, 1.852) * std::pow(p.diameter, 4.87));
  //
  //     // Add nonlinear term: -K * q * |q|^0.852
  //     // Minotaur approximates |q|^0.852 as pow(q^2,0.426)
  //     nlf->addTerm(qij, -K);  // linear placeholder
  //     nlf->setNonlinearOperator(
  //         qij, pow(qij->getSolutionValue() * qij->getSolutionValue(),
  //                  0.426));  // approximate
  //
  //     // Add constraint
  //     p_->newConstraint(nlf, 0.0, 0.0, "headloss_" + std::to_string(a.startNode) + "_" + std::to_string(a.endNode) + "_" + std::to_string(p.id));
  //   }
  // }

  std::cout << me_ << "All constraints successfully added to the problem.\n";
}

void Wdn::solve()
{
  e_->load(p_);
  e_->solve();

  std::cout << "Status: " << e_->getStatusString() << std::endl;
  std::cout << "Objective value = " << e_->getSolutionValue() << std::endl;

  for (HVarConstIter hvar_iter = hvar_.cbegin(); hvar_iter != hvar_.cend(); ++hvar_iter)
    std::cout << hvar_iter->second->getName() << " = " << e_->getSolution()->getPrimal()[hvar_iter->second->getIndex()] << std::endl;

  for (QVarConstIter qvar_iter = qvar_.cbegin(); qvar_iter != qvar_.cend(); ++qvar_iter){
    std::cout<< qvar_iter->second->getName() << " = " << e_->getSolution()->getPrimal()[qvar_iter->second->getIndex()] << std::endl; 
  }
  for (LVarConstIter lvar_iter = lvar_.cbegin(); lvar_iter != lvar_.cend(); ++lvar_iter) {
    std::cout << lvar_iter->second->getName() << " = " << e_->getSolution()->getPrimal()[lvar_iter->second->getIndex()] << std::endl;
  }
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

void Wdn::showHelp() const
{
  std::cout << "Global optimization for water distribution network design problems\n"
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
        << "Global optimization for water distribution network design problems\n";
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
