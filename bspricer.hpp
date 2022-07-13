#ifndef BSPRICER_HPP
#define BSPRICER_HPP

#include <vector>
#include <cmath>


void print(const std::vector<double>& v);


void test_TridiagonalMatrix();

class TridiagonalMatrix {

	// tridiagonal square matrix with a lower diagonal, b diagonal, c upper diagonal
	
public:
	TridiagonalMatrix(int m) : n(m) {}
	
	std::vector<double>& get_a() { return a; }
	std::vector<double>& get_b() { return b; }
	std::vector<double>& get_c() { return c; }
	
	// return x s.t. A * x = rhs
	void solve(
		const std::vector<double>& rhs,
		std::vector<double>& x);

	// returh rhs = A * v
	void product(
		const std::vector<double>& v,
		std::vector<double>& rhs) const;
	
private:
	void lu_decomposte();
	void lu_solve(
		const std::vector<double>& rhs,
		std::vector<double>& x) const;

	int n;
	std::vector<double> a, b, c, l, u, d;
};


class BSVanillaPricer {
	
	// finite difference pricer as per bsv
	
public:
	double calculate(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta,		// 1 for call, -1 for put
		int accuracy);
		
	void test();
	
protected:
	// space grids
	void set_space_grids(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		std::vector<double>& s,
		size_t& idx);

	void set_space_range(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double& smin,
		double& smax);

	void set_concentrated_grids(
		double smin,
		double smax,
		double sfixed,
		std::vector<double>& s);

	// adjust closed pt so that x0 is on grid
	size_t adjust_grids(
		double x0, 
		std::vector<double>& s) const;
	
	// finite difference schemes
	void set_aux(
		const std::vector<double>& s,
		std::vector<std::vector<double> >& aux);
		
	void set_matrix_rhs_fim(
		const std::vector<double>& s,
		double r,
		double q,
		double sigma,
		double dt,
		TridiagonalMatrix& lhs_mat,
		TridiagonalMatrix& rhs_mat);
		
	void set_matrix_rhs_cn(
		const std::vector<double>& s,
		double r,
		double q,
		double sigma,
		double dt,
		TridiagonalMatrix& lhs_mat,
		TridiagonalMatrix& rhs_mat);
	
	// bs formula for sanity checks
	double bspricer(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta) const;
	
	// N(x)
	double normalCDF(double x) const
	{
		return erfc(-x / std::sqrt(2)) / 2;
	}
		
};

#endif