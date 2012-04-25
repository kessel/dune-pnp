#include"example02_operator.hh"
#include<dune/pdelab/localoperator/idefault.hh>

/** a local operator for solving the equation
 *
 *  \partial_t u - \Delta u + a*u = f   in \Omega
 *                              u = g   on \Gamma_D\subseteq\partial\Omega
 *             - \nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *                              u = u_0 at t=t_0
 *
 * (spatial part!) with conforming finite elements on all types of grids in any dimension
 *
 * \tparam B a function indicating the type of boundary condition
 */
template<class B>
class Example03LocalOperator :
  public Example02LocalOperator<B>,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double> // default methods
{
  B& b;
public:
  Example03LocalOperator (B& b_, unsigned int intorder_=2)
    : Example02LocalOperator<B>(b_,intorder_), b(b_) {}
  void preStep (double time, double dt, int stages) {
    b.setTime(time); // enable change of boundary condition type
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep(time,dt,stages);
  }
};
