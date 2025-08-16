#ifndef WDN_H
#define WDN_H

#include "Problem.h"
#include "Objective.h"
#include "Engine.h"
#include "Types.h"
#include <string>
#include <map>
#include <vector>

namespace Minotaur
{
class Wdn {
public:
    Wdn(EnvPtr env);
    ~Wdn();

    /// Starting setup for glob
    void doSetup();

    /// show help messages
    void showHelp() const;

    /// Display information
    int showInfo();


    struct Node {
        int id;
        double elev;
        double demand;
        double pressureMin;
        double pressureMax;
    };
    
    struct Arc {
        int i, j;
        double length;
        double vmax;
        int pipeId; // associated pipe type
    };
    
    struct Pipe {
        int id;
        double diameter;
        double cost;
        double roughness;
    };

    void loadData(const std::string& fileName);
    void printData() const;
    void buildModel();
    void solve();

private:
    const static std::string me_;
    EnvPtr env_;
    ProblemPtr p_;
    EnginePtr e_;

    std::vector<Node> nodes_;
    std::vector<Arc> arcs_;
    std::vector<Pipe> pipes_;
    std::vector<int> sources_;   // <-- multiple sources

    // Decision variables
    std::map<int, VariablePtr> hvar_;                  // h[i]: head at node i
    std::map<std::pair<int,int>, VariablePtr> qvar_;   // q[i,j]: flow on arc (i,j)
    std::map<std::tuple<int,int,int>, VariablePtr> lvar_;  // l[i,j,k]: length of pipe diameter k in arc (i,j)

    void addObjective();
    void addConstraints();
    void setInitialOptions_();
};
}
#endif
