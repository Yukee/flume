#ifndef EULERSOLVER_H
#define EULERSOLVER_H

#include "timeSolver.h"

class EulerSolver : public timeSolver
{
    public:
        EulerSolver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions);

        virtual void get_solution(std::string name, double dt);
};

#endif // EULERSOLVER_H
