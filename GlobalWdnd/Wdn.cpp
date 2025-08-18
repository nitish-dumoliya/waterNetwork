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
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>

using namespace Minotaur;
using namespace std;

const std::string Wdn::me_ = "Wdnd-Global: ";

Wdn::Wdn(EnvPtr env)
  : env_(env)
{
  p_ = (ProblemPtr) new Problem(env_);
  // e_ = ((OsiLPEnginePtr) new OsiLPEngine(env_));
  e_ = ((IpoptEnginePtr) new IpoptEngine(env_));
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

void Wdn::loadData(const string &fileName)
{
  ifstream in(fileName);
  if (!in) {
    cerr << "Error: cannot open file " << fileName << endl;
    return;
  }

  string line;
  while (getline(in, line)) {
    if (line.find("param: nodes:") != string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        istringstream iss(line);
        Node n;
        iss >> n.id >> n.elev >> n.demand >> n.pressureMin >> n.pressureMax;
        if (iss)
          nodes_[n.id] = n;  // insert into map
      }
    } else if (line.find("param: arcs:") != string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        istringstream iss(line);
        Arc a;
        iss >> a.startNode >> a.endNode >> a.length >> a.vmax;
        if (iss)
          arcs_[{a.startNode, a.endNode}] = a;  // insert into map
      }
    } else if (line.find("param: pipes:") != string::npos) {
      while (getline(in, line) && !line.empty() && line[0] != '#') {
        istringstream iss(line);
        Pipe p;
        iss >> p.id >> p.diameter >> p.cost >> p.roughness;
        if (iss)
          pipes_[p.id] = p;  // insert into map
      }
    } else if (line.find("set Source:=") != string::npos) {
      // Example: set Source := 1 4 7 ;
      istringstream iss(line);
      string tmp;
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
  for (const auto &[id, n] : nodes_) {
    std::cout << "  id=" << n.id << " elev=" << n.elev
              << " demand=" << n.demand << " pmin=" << n.pressureMin
              << " pmax=" << n.pressureMax << "\n";
  }

  std::cout << "\nArcs:\n";
  for (const auto &[key, a] : arcs_) {
    std::cout << "  " << a.startNode << " -> " << a.endNode
              << " length=" << a.length << " vmax=" << a.vmax << "\n";
    // << " pipeId=" << a.pipeId << "\n";
  }

  std::cout << "\nPipes:\n";
  for (const auto &[id, p] : pipes_) {
    std::cout << "  id=" << p.id << " d=" << p.diameter << " C=" << p.cost
              << " R=" << p.roughness << "\n";
  }

  std::cout << "\nSource nodes:";
  for (const auto &s : sources_)
    std::cout << " " << s;
  std::cout << "\n";
}

double Wdn::calculateQmax()
{
  double Qmax = 0.0;
  for (const auto &[id, n] : nodes_) {
    if (sources_.find(id) == sources_.end())
      Qmax += n.demand;
  }
  cout << me_ << "Calculated Qmax = " << Qmax << endl;
  return Qmax;
}

void Wdn::buildModel()
{
  OptionDBPtr options = env_->getOptions();
  double Qmax = calculateQmax();

  // Decision variables h[i]
  for (std::map<int, Node>::iterator it = nodes_.begin(); it != nodes_.end();
       ++it) {
   int id = it->first;
    Node &n = it->second;
    std::string vname = "h_" + std::to_string(id);
    VariablePtr h = p_->newVariable(
        n.elev + n.pressureMin, n.elev + n.pressureMax, Continuous, vname);
    hvar_[id] = h;
  }
  std::cout << me_ << "Decision variable h[i] created for " << nodes_.size()
            << " nodes.\n";

  // Decision variables q[i,j]
  for (std::map<std::pair<int, int>, Arc>::iterator it = arcs_.begin();
       it != arcs_.end(); ++it) {
    std::pair<int, int> key = it->first;
    Arc &a = it->second;
    VariablePtr q = p_->newVariable(
        -Qmax, Qmax, Continuous,
        "q_" + std::to_string(a.startNode) + "_" + std::to_string(a.endNode));
    qvar_[key] = q;
  }
  std::cout << me_ << "Decision variable q[i,j] created for " << arcs_.size()
            << " arcs.\n";

  // Decision variables l[i,j,k]
  for (std::map<std::pair<int, int>, Arc>::iterator it = arcs_.begin();
       it != arcs_.end(); ++it) {
    Arc &a = it->second;

    for (std::map<int, Pipe>::iterator pit = pipes_.begin();
         pit != pipes_.end(); ++pit) {
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

  for (const auto &[key, a] : arcs_) {
    for (const auto &[pid, p] : pipes_) {
      VariablePtr l = lvar_[{a.startNode, a.endNode, pid}];
      lf->addTerm(l, p.cost * a.length);
    }
  }

  FunctionPtr f = (FunctionPtr) new Function(lf);
  obj_ = p_->newObjective(f, 0.0, Minimize);
  cout << me_ << "Objective function successfully added to the problem.\n";
}

void Wdn::addConstraints()
{
  /*
   * Constraints:
   *
   * 1. Flow balance at nodes:
   *      ∑_{j: (i,j) ∈ arcs} q[i,j] - ∑_{j: (j,i) ∈ arcs} q[j,i] =
   * demand[i]
   *
   *    Ensures conservation of flow: inflow - outflow = demand.
   *
   * 2. Head-loss (Hazen-Williams simplified):
   *      h[i] - h[j] = K * q[i,j]*|q[i,j]|^0.852
   *
   *    with K = 10.67 * L / (C^1.852 * d^4.87)
   *    Here we approximate nonlinearities with NLP terms.
   */

  // 1. Flow balance
  for (const auto &[id, n] : nodes_) {
    if (sources_.find(id) != sources_.end())
      continue;
    LinearFunctionPtr lf = new LinearFunction();
    for (const auto &[key, a] : arcs_) {
      if (a.startNode == id)
        lf->addTerm(qvar_[key], -1.0);
      if (a.endNode == id)
        lf->addTerm(qvar_[key], 1.0);
    }
    FunctionPtr f = (FunctionPtr) new Function(lf);
    p_->newConstraint(f, n.demand, n.demand, "flow_balance_" + to_string(id));
  }

  // 2. Pipe lengths = arc length
  for (const auto &[key, a] : arcs_) {
    LinearFunctionPtr lf = new LinearFunction();
    for (const auto &[pid, p] : pipes_) {
      lf->addTerm(lvar_[{a.startNode, a.endNode, pid}], 1.0);
    }
    FunctionPtr f = (FunctionPtr) new Function(lf);
    p_->newConstraint(
        f, a.length, a.length,
        "con3_" + to_string(a.startNode) + "_" + to_string(a.endNode));
  }

  // 3. Fix head at source nodes
  for (const auto &s : sources_) {
    const Node &srcNode = nodes_.at(s);
    LinearFunctionPtr lf = new LinearFunction();
    lf->addTerm(hvar_[s], 1.0);
    FunctionPtr f = (FunctionPtr) new Function(lf);
    p_->newConstraint(f, srcNode.elev + srcNode.pressureMin,
                      srcNode.elev + srcNode.pressureMin,
                      "con4_headFix_" + to_string(s));
  }

  // 4. Head-loss constraints(Hazen-Williams formula)
  // for (std::map<std::pair<int, int>, Arc>::iterator ait = arcs_.begin();
  //      ait != arcs_.end(); ++ait) {
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
  //   for (std::map<int, Pipe>::iterator pit = pipes_.begin();
  //        pit != pipes_.end(); ++pit) {
  //     Pipe &p = pit->second;
  //
  //     // Get length variable
  //     VariablePtr lijk = lvar_[{a.startNode, a.endNode, p.id}];
  //
  //     // Compute Hazen-Williams coefficient
  //     double K = 10.67 * a.length /
  //                (std::pow(p.roughness, 1.852) *
  //                std::pow(p.diameter, 4.87));
  //
  //
  //     // Add nonlinear term: -K * q * |q|^0.852
  //     // Minotaur approximates |q|^0.852 as pow(q^2,0.426)
  //     nlf->addTerm(qij, -K);  // linear placeholder
  //     nlf->setNonlinearOperator(
  //         qij, pow(qij->getSolutionValue() * qij->getSolutionValue(),
  //                  0.426));  // approximate
  //
  //     // Add constraint
  //     p_->newConstraint(nlf, 0.0, 0.0,
  //                       "headloss_" + std::to_string(a.startNode) + "_" +
  //                           std::to_string(a.endNode) + "_" +
  //                           std::to_string(p.id));
  //   }
  // }

  // === Head-loss constraints (Hazen-Williams) ===
  
  // for (std::map<std::pair<int, int>, Arc>::iterator ait = arcs_.begin();
  //      ait != arcs_.end(); ++ait) {
  //
  //   int i = ait->second.startNode;
  //   int j = ait->second.endNode;
  //
  //   CGraph *nlf = new CGraph();
  //   CNode *c_h_i, *c_h_j, *c_q, *c_absq, *c_exp, *c_pipeSum, *c_tmp;
  //
  //   // Node heads
  //   c_h_i = nlf->newNode(hvar_[i]);
  //   c_h_j = nlf->newNode(hvar_[j]);
  //
  //   // Flow q[i,j]
  //   c_q = nlf->newNode(qvar_[{i, j}]);
  //
  //   // Absolute value: |q|
  //   // c_absq = nlf->newNode(OpAbs, c_q);
  //   CNode *cq_sq = nlf->newNode(OpMult, c_q, c_q);  // q^2
  //   CNode *cq_abs =
  //       nlf->newNode(OpPowK, cq_sq, nlf->newNode(0.5));  // sqrt(q^2)
  //
  //   // |q|^0.852
  //   c_exp = nlf->newNode(0.852);
  //   c_absq = nlf->newNode(OpPowK, cq_abs, c_exp);
  //
  //   // q * |q|^0.852
  //   c_q = nlf->newNode(OpMult, c_q, c_absq);
  //
  //   // Build sum over pipes
  //   c_pipeSum = nlf->newNode(0.0);  // initialize
  //
  //   // for (std::vector<Pipe>::iterator pit = pipes_.begin();
  //   // pit != pipes_.end(); ++pit)
  //   for (std::map<int, Pipe>::iterator pit = pipes_.begin();
  //        pit != pipes_.end(); ++pit) {
  //
  //     int k = pit->first;
  //     const Pipe &pipe = pit->second;
  //     // l[i,j,k] (decision variable)
  //     CNode *c_l = nlf->newNode(lvar_[{i, j, k}]);
  //
  //     // coefficient = 10.67 / (R^1.852 * D^4.87)
  //     double coeff = 10.67 / (std::pow(pipe.roughness, 1.852) *
  //                             std::pow(pipe.diameter, 4.87));
  //
  //     CNode *c_coeff = nlf->newNode(coeff);
  //
  //     // term = l[i,j,k] * coeff
  //     c_tmp = nlf->newNode(OpMult, c_l, c_coeff);
  //
  //     // add to sum
  //     c_pipeSum = nlf->newNode(OpPlus, c_pipeSum, c_tmp);
  //   }
  //
  //   // Multiply flow term by pipe sum
  //   CNode *c_loss = nlf->newNode(OpMult, c_q, c_pipeSum);
  //
  //   // Constraint: (h_i - h_j) - loss = 0
  //   CNode *c_diff = nlf->newNode(OpMinus, c_h_i, c_h_j);
  //   CNode *c_final = nlf->newNode(OpMinus, c_diff, c_loss);
  //
  //   nlf->setOut(c_final);
  //   nlf->finalize();
  //
  //   FunctionPtr f = (FunctionPtr) new Function(nlf);
  //   p_->newConstraint(f, 0.0, 0.0,
  //                     "HL_" + std::to_string(i) + "_" + std::to_string(j));
  // }
  cout << me_ << "All constraints successfully added to the problem.\n";
}

void Wdn::solve()
{
  e_->load(p_);
  e_->solve();

  cout << "Status: " << e_->getStatusString() << endl;
  cout << "Objective value = " << e_->getSolutionValue() << endl;

  // for (const auto &[id, var] : hvar_)
  //   cout << "Head at node " << id << " = "
  //        << e_->getSolution()->getPrimal()[var->getIndex()] << endl;
  //
  // for (const auto &[key, var] : qvar_)
  //   cout << "Flow on arc (" << key.first << "," << key.second
  //        << ") = " << e_->getSolution()->getPrimal()[var->getIndex()] <<
  //        endl;
  //
  // for (const auto &[key, var] : lvar_) {
  //   cout << "l(" << get<0>(key) << "," << get<1>(key) << "," << get<2>(key)
  //        << ") = " << e_->getSolution()->getPrimal()[var->getIndex()] <<
  //        endl;
  // }
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
  cout << "Global optimization for water distribution network design "
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
    options->write(cout);
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
        << me_ << "Minotaur version " << env_->getVersion() << endl
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
      << me_ << "Minotaur version " << env_->getVersion() << endl
      << me_
      << "Global optimization for water distribution network design "
         "problems\n";
  return 0;
}
