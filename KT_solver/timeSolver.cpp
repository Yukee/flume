#include "timeSolver.h"
#include <iostream>
#include <fstream>

using namespace std;

timeSolver::timeSolver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions)
: m_deltaT(deltaT), m_T(T), m_spatialSolver(spatialSolver)
{
    m_un = initial_conditions; //Why can't I use the initialization list?
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
	infos.open( ("Results/" + path + "_infos.tsv").c_str(), ios::out);

	infos << "Writing results in files named " << path + "_i.tsv" << endl;

	infos << "timestep\t" << m_dt << endl;

	infos << "endtime\t" << m_T << endl;

	int n = m_spatialSolver->get_space_dimensions();
	infos << "spacedimensions\t" << n << endl;

	infos << "solveddimensions\t" << m_spatialSolver->get_solved_dimensions() << endl;

	Vector< Vector<double> > dom = m_spatialSolver->get_domain_bounds();
	for(int d=0;d<n;d++) 
	{
		infos << "domain\t" << d << "\t" << dom[0][d] << "\t" << dom[1][d] << endl;
	}

	Vector<double> dx = m_spatialSolver->get_deltaX();
	for(int d=0;d<n;d++) 
	{
		infos << "cellsdimension\t" << d << "\t" << dx[d] << endl;
	}
}
