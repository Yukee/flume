#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <vector>
#include <string>
#include <fstream>
#include "FD1Solver.h"
#include "ScalarField.h"

class timeSolver
{
    public:
        timeSolver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions);

	// records results in a file every every dt seconds (in simulation time)
        virtual void get_solution(std::string path, double dt);

	// records numbFiles results files regularly in the give integration time m_T
	virtual void get_solution(std::string path, int numbFiles);

	// print infos about the simulation in a file	
	void print_infos(std::string path);

    protected:
        double m_deltaT;
        double m_T;
	double m_dt;
        int m_ntSteps;
        FD1Solver *m_spatialSolver;
        VectorField m_un;
	
};

#endif // TIMESOLVER_H
