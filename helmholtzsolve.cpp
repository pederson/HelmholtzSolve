#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/SparseVector.hpp"
#include "LinearAlgebra/SparseMatrix.hpp"
#include "LinearAlgebra/Preconditioner.hpp"

using namespace std;

// g++ -std=c++11 fishsolve.cpp -o fishsolve

#include <typeinfo>
#include "LinearAlgebra/LinearSolvers.hpp"

int main(int argc, char * argv[]){

	unsigned int ny=100, nx = 100;
	double xmin = 0, xmax = 1;
	double ymin = 0, ymax = 1;

	double dx = (xmax-xmin)/double(nx-1);
	double dy = (ymax-ymin)/double(ny-1);

	double k = dx/2;
	double eps0 = 8.854e-12;
	double mu0 = 4*3.14159e-7;

	// construct the operator matrix
	SparseMatrix Op(nx*ny, nx*ny);
	unsigned int c, up, dn, lt, rt;
	for (auto i=0; i<nx; i++){
		for (auto j=0; j<ny; j++){
			c = nx*j + i;
			up = nx*(j+1) + i;
			dn = nx*(j-1) + i;
			lt = nx*j + i-1;
			rt = nx*j + i+1;

			Op.set(c, c, 4 - 0.75);

			if (j<ny-1)	Op.set(c, up, -1);
			if (j>0) 	Op.set(c, dn, -1);
			if (i>0)	Op.set(c, lt, -1);
			if (i<nx-1) Op.set(c, rt, -1);
		}
	}

	// right hand side forcing
	Vector rho(nx*ny);
	rho(nx*(ny*0.5) + nx*0.5) = 1*dx*dx;
	/*
	// quadrupole
	rho(nx*(ny*0.25)+nx*0.25) = -1*dx*dx;
	rho(nx*(ny*0.25)+nx*0.75) = 1*dx*dx;
	rho(nx*(ny*0.75)+nx*0.25) = 1*dx*dx;
	rho(nx*(ny*0.75)+nx*0.75) = -1*dx*dx;
	*/



	cout << "****************** PRECONDITIONED SOLVERS ******************" << endl;
	Vector x(nx*ny);
	cout << "nnz: " << Op.nnz() << "/" << Op.rows()*Op.cols() << " = " << double(Op.nnz())/double(Op.rows()*Op.cols()) << endl;
	unsigned int niters;
	unsigned int nitmax = 5000;
	double tol = 1e-14;

	// cout << "************************** CG - JACOBI PC:" << endl;
	// JacobiPreconditioner jpc(Op);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&jpc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;


	// cout << "************************** CG - GAUSS-SEIDEL PC:" << endl;
	// GSPreconditioner gspc(Op);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&gspc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;


	// cout << "************************** CG - SYMMETRIC GAUSS-SEIDEL PC:" << endl;
	// SGSPreconditioner sgspc(Op);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&sgspc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;


	// cout << "************************** CG - SOR PC:" << endl;
	// SORPreconditioner sorpc(Op, 1.5);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&sorpc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;


	// cout << "************************** CG - SSOR PC:" << endl;
	// SSORPreconditioner ssorpc(Op, 0.8);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&ssorpc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;


	// cout << "************************** CG - INCOMPLETE CHOLESKY PC:" << endl;
	// ICPreconditioner icpc(Op);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&icpc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;

	// cout << "************************** CG - INCOMPLETE LU PC:" << endl;
	// ILUPreconditioner ilupc(Op);
	// x.fill(0);
	// niters = conjugate_gradient(Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&ilupc, Op, rho, x, nitmax);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << (rho - Op*x).norm() << endl;

	cout << "************************** CG - AMG PC:" << endl;
	// AMGPreconditioner amgpc(Op);
	x.fill(0);
	niters = conjugate_gradient(Op, rho, x, nitmax, tol);
	cout << "iterated: " << niters << " times" << endl;
	cout << "cg resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = conjugate_gradient(&amgpc, Op, rho, x, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc cg resid: " << norm_2(rho - Op*x) << endl;


	// cout << "************************** BICGSTAB - INCOMPLETE LU PC:" << endl;
	// x.fill(0);
	// niters = bicgstab(Op, rho, x, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "bicgstab resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = bicgstab(&ilupc, Op, rho, x, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc bicgstab resid: " << norm_2(rho - Op*x) << endl;

	// cout << "************************** GMRES - INCOMPLETE LU PC:" << endl;
	// x.fill(0);
	// niters = gmres_k(Op, rho, x, 20, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "gmres resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = gmres_k(&ilupc, Op, rho, x, 20, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc gmres resid: " << norm_2(rho - Op*x) << endl;

	// cout << "************************** GMRES - AMG PC:" << endl;
	// x.fill(0);
	// niters = gmres_k(Op, rho, x, 20, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "gmres resid: " << (rho - Op*x).norm() << endl;
	// x.fill(0);
	// niters = gmres_k(&amgpc, Op, rho, x, 20, nitmax, tol);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "pc gmres resid: " << norm_2(rho - Op*x) << endl;



	x.dlmwrite("solution.txt");



	// cout << "************************** ALGEBRAIC MULTIGRID: " << endl;
	// Vector offd(nitmax - 1); offd.fill(-1);
	// SparseMatrix spamg = 2*speye(nitmax,nitmax) + spdiag(offd, 1) + spdiag(offd,-1);
	// Vector ps_amg(nitmax); ps_amg.fill(0);
	// Vector psamg_b = randvecn(nitmax);
	// cout << "amg resid before: " << (psamg_b - spamg*ps_amg).norm() << endl;
	// niters = amg(spamg, psamg_b, ps_amg, nitmax, 1.0e-12);
	// cout << "iterated: " << niters << " times" << endl;
	// cout << "amg resid: " << (psamg_b - spamg*ps_amg).norm() << endl;
	// ps_amg.dlmwrite("amg_soln.txt");

	return 0;
}
