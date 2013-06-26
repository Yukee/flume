#include "timeSolver.h"
#include <iostream>
#include <fstream>

using namespace std;

timeSolver::timeSolver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions)
: m_deltaT(deltaT), m_T(T), m_spatialSolver(spatialSolver), m_un(initial_conditions)
{
    m_ntSteps = m_T/m_deltaT;
}

void timeSolver::get_solution(string path, double dt)
{
	m_dt = dt;
	print_infos(path);
}

void timeSolver::get_solution(string path, int numbFiles)
{
	get_solution(path, m_T/numbFiles);
}

void timeSolver::print_infos(string path)
{
	fstream infos;
	infos.open( ("Results/" + path + "_standard_infos.tsv").c_str(), ios::out);

	int n = m_spatialSolver->get_space_dimensions();
	Vector< Vector<double> > dom = m_spatialSolver->get_domain_bounds();
	Vector<int> nb_cells = m_spatialSolver->get_nxSteps();

	infos << "Writing results in files named " << path + "_i.tsv" << endl;

	infos << m_dt << endl << m_T << endl;

	for(int d=0;d<n;d++)
	{
		infos << dom[0][d] << endl << dom[1][d] << endl;
	}
	
	for(int d=0;d<n;d++)
	{
		infos << nb_cells[d] << endl;
	}
}
