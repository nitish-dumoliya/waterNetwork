//
//    Minotaur -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2025 The Minotaur Team.
//

/**
 * \file WdnMain.cpp
 * \brief The main function for solving instances with wdnd
 * \author Nitish Kumar Dumoliya, IIT Bombay
 */

#include "MinotaurConfig.h"
#include "Environment.h"
#include "Wdn.h"
#include "Logger.h"
#include "Option.h"
#include "Problem.h"
#include "Types.h"

#include <iostream>

using namespace Minotaur;

int main(int argc, char **argv)
{
  EnvPtr env = (EnvPtr) new Environment();
  ProblemPtr p = (ProblemPtr) new Problem(env) ;
  std::string datfile;
  int err = 0;
  Wdn wdn(env);

  wdn.doSetup();

  // read user-specified options
  env->readOptions(argc, argv);
  // env->getOptions()->findBool("use_native_cgraph")->setValue(false);
  
  if (argc > 1) {
    datfile = argv[1];
  }

  if (0 != wdn.showInfo()) {
    goto CLEANUP;
  }

  wdn.loadData(datfile);
  wdn.printData();
  p = wdn.buildModel(p);
  err = wdn.solve(p);
  if (err) {
      goto CLEANUP;
  }

CLEANUP:
  if (p) {
    delete p;
  }
  if (env) {
    delete env;
  }

  return 0;
}
