#include "Wdn.h"
#include "EngineFactory.h"
#include "Environment.h"
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
#include "Objective.h"
#include "Option.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace Minotaur;
const std::string Wdn::me_ = "Wdnd-Global: ";

using namespace std;

Wdn::Wdn(EnvPtr env) : env_(env)
{
    p_ = (ProblemPtr) new Problem(env_);
    e_ = EngineFactory(env_).getNLPEngine(); // IPOPT
}

Wdn::~Wdn()
{ 
  // if (sol_) {
  //   delete sol_;
  // }
}

void Wdn::doSetup()
{
  setInitialOptions_();
}

void Wdn::loadData(const string& fileName) {
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
                if (iss) nodes_.push_back(n);
            }
        }
        else if (line.find("param: arcs:") != string::npos) {
            while (getline(in, line) && !line.empty() && line[0] != '#') {
                istringstream iss(line);
                Arc a;
                iss >> a.i >> a.j >> a.length >> a.vmax;
                if (iss) arcs_.push_back(a);
            }
        }
        else if (line.find("param: pipes:") != string::npos) {
            while (getline(in, line) && !line.empty() && line[0] != '#') {
                istringstream iss(line);
                Pipe p;
                iss >> p.id >> p.diameter >> p.cost >> p.roughness;
                if (iss) pipes_.push_back(p);
            }
        }
	    else if (line.find("set Source:=") != string::npos) {
            // Example: set Source := 1 4 7 ;
            istringstream iss(line);
            string tmp;
            iss >> tmp; // "set"
            iss >> tmp; // "Source:="
            int s;
            while (iss >> s) {
                sources_.push_back(s);
            }
        }
    }
}

void Wdn::printData() const {
    std::cout << "Nodes:\n";
    for (const auto& n : nodes_) {
        std::cout << "  id=" << n.id << " elev=" << n.elev
                  << " demand=" << n.demand
                  << " pmin=" << n.pressureMin
                  << " pmax=" << n.pressureMax << "\n";
    }

    std::cout << "\nArcs:\n";
    for (const auto& a : arcs_) {
        std::cout << "  " << a.i << " -> " << a.j
                  << " length=" << a.length << " vmax=" << a.vmax << "\n";
    }

    std::cout << "\nPipes:\n";
    for (const auto& p : pipes_) {
        std::cout << "  id=" << p.id << " d=" << p.diameter
                  << " C=" << p.cost
                  << " R=" << p.roughness << "\n";
    }
    std::cout << "\nSource nodes:";
    for (int s : sources_) std::cout << " " << s;
    std::cout << "\n";
}

void Wdn::buildModel()
{
    OptionDBPtr options = env_->getOptions();
    // === Decision variables ===
    for (auto &n : nodes_) {
        Node node = n;
        VariablePtr h = p_->newVariable(node.pressureMin, node.pressureMax,
                                Continuous, "h"+std::to_string(node.id));
        hvar_[node.id] = h;
    }
    std::cout << me_ << "Decision variable h[i] created for "
              << nodes_.size() << " nodes.\n";

    for (auto &a : arcs_) {
        VariablePtr q = p_->newVariable(-a.vmax, a.vmax, Continuous,
                            "q_"+std::to_string(a.i)+"_"+std::to_string(a.j));
        qvar_[{a.i,a.j}] = q;
    }
    std::cout << me_ << "Decision variable q[i,j] created for "
              << arcs_.size() << " arcs.\n";

    for (auto &a : arcs_) {
        for (auto &p : pipes_) {
            std::string vname = "l_" + std::to_string(a.i) + "_" +
                                          std::to_string(a.j) + "_" +
                                          std::to_string(p.id);
            // Binary: either pipe k is selected on arc (i,j) or not
            VariablePtr l = p_->newVariable(0.0, a.length, Continuous, vname);
            lvar_[{a.i, a.j, p.id}] = l;
        }
    }
    std::cout << me_ << "Decision variable l[i,j,k] created for "
              << arcs_.size() << " arcs × "
              << pipes_.size() << " pipes.\n";
    // === Objective Function ===
    addObjective();
    // === Constraints ===
    addConstraints();
    if(options->findBool("display_problem")->getValue() == true) {
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
    * This ensures that the chosen pipe design minimizes the total installation 
    * cost.
    */
    ObjectivePtr obj_;
    VariablePtr l=0;
    FunctionPtr f=0;
    LinearFunctionPtr lf = new LinearFunction();
    LinearFunctionPtr lf2 = new LinearFunction();

    for (auto &a : arcs_) {
        for (auto &p : pipes_) {
            // Get variable l[i,j,k]
            l = lvar_[{a.i, a.j, p.id}];
            // Add term: cost[k] * length(i,j) * l[i,j,k]
	        lf->addTerm(l, p.cost * a.length);
            lf2->add(lf); 
        }
    }

    f = (FunctionPtr) new Function(lf2);
    obj_ = p_->newObjective(f,0.0, Minimize);
    std::cout << me_ <<  "Objective function successfully added to the problem." 
                     << std::endl;
    // obj_->write(std::cout); 
}


void Wdn::addConstraints()
{
    /*
     * Constraints:
     *
     * 1. Flow balance at nodes:
     *      ∑_{j: (i,j) ∈ arcs} q[i,j] - ∑_{j: (j,i) ∈ arcs} q[j,i] = demand[i]
     *
     *    Ensures conservation of flow: inflow - outflow = demand.
     *
     * 2. Head-loss (Hazen-Williams simplified):
     *      h[i] - h[j] = K * q[i,j]*|q[i,j]|^0.852
     *
     *    with K = 10.67 * L / (C^1.852 * d^4.87)
     *    Here we approximate nonlinearities with NLP terms.
     */

    // === 1. Flow balance at each node ===
    for (auto &n : nodes_) {
        LinearFunctionPtr lf = new LinearFunction();
        Node node = n;

        for (auto &a : arcs_) {
            if (a.i == node.id) {
                lf->addTerm(qvar_[{a.i, a.j}], -1.0); // outflow
            }
            if (a.j == node.id) {
                lf->addTerm(qvar_[{a.i, a.j}], +1.0); // inflow
            }
        }

        // equality: inflow - outflow = demand
        p_->newConstraint(lf, node.demand, node.demand, 
                          "flow_balance_" + std::to_string(n.id));
    }

    // === 2. Head-loss constraints ===
    // for (auto &a : arcs_) {
    //     Node ni = nodes_[a.i];
    //     Node nj = nodes_[a.j];
    //
    //     // All candidate pipes for this arc (if multiple diameters are possible)
    //     for (auto &p : pipes_) {
    //         std::array<int,3> key = {a.i, a.j, p.id};
    //         if (lvar_.find(key) == lvar_.end()) continue; // skip if pipe var not defined
    //
    //         VariablePtr q  = qvar_[{a.i, a.j}];
    //         VariablePtr hi = hvar_[ni.id];
    //         VariablePtr hj = hvar_[nj.id];
    //
    //         double C = p.roughness;
    //         double d = p.diameter;
    //
    //         // Hazen-Williams coefficient
    //         double K = 10.67 * a.length / (std::pow(C, 1.852) * std::pow(d, 4.87));
    //
    //         // Head loss expression: hi - hj = K * q|q|^0.852
    //         // Minotaur allows nonlinear constraints using Function objects
    //         NonlinearFunctionPtr nlf = new NonlinearFunction();
    //
    //         // left side: hi - hj
    //         nlf->addTerm(hi, 1.0);
    //         nlf->addTerm(hj, -1.0);
    //
    //         // right side: -K * |q|^1.852
    //         // NOTE: this requires adding a nonlinear expression, not just QuadraticFunction
    //         // For demonstration, approximate q*|q|^0.852
    //         nlf->addTerm(q, -K);  // simplified linearization
    //         // TODO: replace with exact nonlinear operator if solver supports pow(q, 1.852)
    //
    //         p_->newConstraint(nlf, 0.0, 0.0,
    //                           "headloss_" + std::to_string(a.i) + "_" + std::to_string(a.j));
    //     }
    // }

    std::cout << me_ << "All constraints successfully added to the problem." << std::endl;
}


// void Wdn::addConstraints()
// {
//     // === Flow balance at nodes ===
//     for (auto &n : nodes_) {
//         LinearFunctionPtr lf = new LinearFunction();
//
//         for (auto &a : arcs_) {
//             if (a.i == n.first) lf->addTerm(qvar_[{a.i,a.j}], -1.0);
//             if (a.j == n.first) lf->addTerm(qvar_[{a.i,a.j}], +1.0);
//         }
//
//         Node node = n.second;
//         p_->newConstraint(lf, node.demand, node.demand, "balance_"+std::to_string(n.first));
//     }
//
//     // === Head-loss constraints (Hazen-Williams equation) ===
//     for (auto &a : arcs_) {
//         Node ni = nodes_[a.i];
//         Node nj = nodes_[a.j];
//         Pipe pipe = pipes_[a.pipeId];
//
//         double C = pipe.roughness;
//         double d = pipe.diameter;
//
//         VariablePtr hi = hvar_[ni.id];
//         VariablePtr hj = hvar_[nj.id];
//         VariablePtr q  = qvar_[{a.i,a.j}];
//
//         // Head loss Δh = K * q|q|^0.852  (Hazen-Williams equation)
//         double K = 10.67 * a.length / (std::pow(C,1.852) * std::pow(d,4.87));
//
//         QuadraticFunctionPtr qf = new QuadraticFunction();
//         qf->addTerm(q, q, K); // NOTE: real Hazen-Williams is nonlinear, here simplified
//
//         LinearFunctionPtr lf = new LinearFunction();
//         lf->addTerm(hi, 1.0);
//         lf->addTerm(hj, -1.0);
//
//         p_->newConstraint(lf, -INFINITY, qf, "headloss_"+std::to_string(a.i)+"_"+std::to_string(a.j));
//     }
// }

// void Wdn::solve()
// {
//     e_->load(p_);
//     e_->solve();
//
//     std::cout << "Status: " << e_->getStatusString() << std::endl;
//     std::cout << "Objective value = " << e_->getSolutionValue() << std::endl;
//
//     for (auto &n : hvar_) {
//         std::cout << "Head at node " << n.first << " = "
//                   << e_->getSolution()->getPrimal()[n.second->getIndex()] << std::endl;
//     }
//
//     for (auto &a : qvar_) {
//         std::cout << "Flow on arc (" << a.first.first << "," << a.first.second << ") = "
//                   << e_->getSolution()->getPrimal()[a.second->getIndex()] << std::endl;
//     }
// }

void Wdn::setInitialOptions_()
{
  OptionDBPtr options = env_->getOptions();
  options->findString("interface_type")->setValue("ampl");
  options->findBool("presolve")->setValue(true);
  options->findBool("nl_presolve")->setValue(true);
  options->findBool("lin_presolve")->setValue(true);
  options->findBool("msheur")->setValue(true);
  // options->findString("brancher")->setValue("hybrid");
  options->findString("nlp_engine")->setValue("ipopt");
  // options->findBool("cgtoqf")->setValue(true);
  // options->findBool("simplex_cut")->setValue(true);
}


void Wdn::showHelp() const
{
  std::cout << "Global optimization for water distribution network design problems" << std::endl
            << "Usage:" << std::endl
            << "To show version: wdnd -v (or --display_version yes)"
            << std::endl
            << "To show all options: wdnd -= (or --display_options yes)"
            << std::endl
            << "To solve an instance: wdnd --option1 [value] "
            << "--option2 [value] ... "
            << " .dat-file" << std::endl;
}


int Wdn::showInfo()
{
  OptionDBPtr options = env_->getOptions();

  if(options->findBool("display_options")->getValue() ||
     options->findFlag("=")->getValue()) {
    options->write(std::cout);
    return 1;
  }

  if(options->findBool("display_help")->getValue() ||
     options->findFlag("?")->getValue()) {
    showHelp();
    return 1;
  }

  if(options->findBool("display_version")->getValue() ||
     options->findFlag("v")->getValue()) {
    env_->getLogger()->msgStream(LogNone)
        << me_ << "Minotaur version " << env_->getVersion() << std::endl
        << me_ << "Global optimization for water distribution network design problems" << std::endl;
    return 1;
  }

  if(options->findString("problem_file")->getValue() == "") {
    showHelp();
    return 1;
  }

  env_->getLogger()->msgStream(LogInfo)
      << me_ << "Minotaur version " << env_->getVersion() << std::endl
      << me_ << "Global optimization for water distribution network design problems" << std::endl;
  return 0;
}
