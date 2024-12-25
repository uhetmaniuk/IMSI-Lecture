#pragma once

namespace msfem {

/// \brief Class to describe the problem to solve
class ProblemDescriptor {

  protected:

    void ReadFile(const char *fName);

  public:

    ProblemDescriptor(const char *fName);

    /// Destructor
    ~ProblemDescriptor() = default;

    /// Solves the problem
    void Solve();

};

} // namespace msfem

