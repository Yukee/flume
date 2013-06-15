#ifndef RK3SOLVER_H
#define RK3SOLVER_H

#include "timeSolver.h"

class RK3Solver : public timeSolver
{
    public:
        RK3Solver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions);

        virtual void get_solution(std::string name, double dt);
};

#endif
