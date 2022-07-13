#ifndef  BS_BARRIER_PRICER_HPP
#define  BS_BARRIER_PRICER_HPP

#include <vector>
#include "bspricer.hpp"


class FDGridHelper {

	// helper to set up space grids with concentrations as lvp
	
public:
	static void 
	set_concentrated_grids(
		double smin,
		double smax,
		const std::vector<double>& critical_points,
		std::vector<double>& vec_space_grid
		);

	static void test();
	
private:		
	static double
	global_jacobian(
        const double y,
        const double alpha,
        const std::vector<double>& critical_points
        );
        
	static void
	ruger_kutta(
        const double x_start,
        const double y_start,
        const double xn,
        const double A,
        const double alpha,
        const std::vector<double>& critical_points,
        std::vector<double>& result
        );
        
	static int
	upper_bound_for_ruger_kutta(
        const double smin,
        const double smax,
        const double alpha,
        const std::vector<double>& critical_points,
        std::vector<double>& vec_space_grid,
        double& high
        );

};


class BSBarrierPricer : public BSVanillaPricer {

	// finite difference pricer as per lvp
	// CN + FIM

public:
	double calculate(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta,		// 1 for call, -1 for put
		double L,
		double U,
		int accuracy);
		
	void test();

protected:
	void cn_step(
		double r,
		double q,
		double sigma,
		double L,
		double U,
		double dt,
		const std::vector<double>& s,
		std::vector<double>& vnext,
		std::vector<double>& v,
		std::vector<double>& rhs,
		bool is_monitored);
		
	void fim_step(
		double r,
		double q,
		double sigma,
		double L,
		double U,
		double dt,
		const std::vector<double>& s,
		std::vector<double>& vnext,
		std::vector<double>& v,
		std::vector<double>& rhs,
		bool is_monitored);
		
	// apply barriers right after FD scheme
	void update_matrix_rhs(
		double L,
		double U,
		const std::vector<double>& s,
		double dt,
		TridiagonalMatrix& lhs_mat,
		TridiagonalMatrix& rhs_mat);
	
	void apply_barriers(
		double L,
		double U,
		const std::vector<double>& s,
		std::vector<double>& v);
		
	// space grids
	void set_space_grids(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double L,
		double U,
		std::vector<double>& s,
		size_t& idx);

};


class BSWindowBarrierPricer : public BSBarrierPricer {

	// finite difference pricer as per lvp

public:
	double calculate(
		double S0,
		double K,
		double r,
		double q,
		double sigma,
		double T,
		double theta,
		double L,
		double U,
		double T1,
		double T2,
		int accuracy);
		
	void test();

protected:
	void set_time_grids(
		double T1,
		double T2,
		double dt,
		std::vector<double>& ts,
		size_t& t1_idx,
		size_t& t2_idx) const;

	bool is_monitored(
		double T1, 
		double T2, 
		double t);

	void apply_barriers(
		double L,
		double U,
		const std::vector<double>& s,
		double T1,
		double T2,
		double t,
		std::vector<double>& v);
		
};
		
#endif