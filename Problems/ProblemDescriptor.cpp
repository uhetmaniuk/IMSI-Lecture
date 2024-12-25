#include "Problems/ProblemDescriptor.h"

namespace msfem {

ProblemDescriptor::ProblemDescriptor(const char *fName) {

  ReadFile(fName);

}

void ProblemDescriptor::Solve() {

//  return (myProblem) ? myProblem->Solve() : -1; 

}

void ProblemDescriptor::ReadFile(const char *fName) {

  ifstream fin(fName);
  if (!fin) {
    std::cerr << "\n\n !! Input data file " << fName << " opening failed !! \n\n\n";
    exit(-1);
  }

}

} // namespace msfem
