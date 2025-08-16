#include "Environment.h"
#include "Wdn.h"
#include "Logger.h"
#include "Option.h"
#include "Problem.h"
#include "Types.h"
#include "MinotaurConfig.h"

#include <iostream>

using namespace Minotaur;

int main(int argc, char **argv)
{
  EnvPtr env = (EnvPtr) new Environment();
  ProblemPtr p_ = 0;
  std::string datfile;
  Wdn wdn(env);

  wdn.doSetup();

  // read user-specified options
  env->readOptions(argc, argv);
  env->getOptions()->findBool("use_native_cgraph")->setValue(true);

  if (argc > 1) {
    datfile = argv[1];
  }

  if (0 != wdn.showInfo()) {
    goto CLEANUP;
  }

  wdn.loadData(datfile);
  wdn.printData();
  wdn.buildModel();
  wdn.solve();
  goto CLEANUP;

CLEANUP:
  if (p_) {
    delete p_;
  }
  if (env) {
    delete env;
  }

  return 0;
}
