#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

/** a local operator for solving the equation
 *
 *   - \Delta u + a*u = f   in \Omega
 *                  u = g   on \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = j   on \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 * with conforming finite elements on all types of grids in any dimension
 *
 * \tparam B a function indicating the type of boundary condition
 */


#define PI 3.1415

template<class PbB, class F, class S, class FluxContainer, class CpDgf, class CmDgf>
class PoissonOperator : 
  public Dune::PDELab::NumericalJacobianApplyVolume       <PoissonOperator<PbB, F, S, FluxContainer, CpDgf, CmDgf> >,
  public Dune::PDELab::NumericalJacobianVolume            <PoissonOperator<PbB, F, S, FluxContainer, CpDgf, CmDgf> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary     <PoissonOperator<PbB, F, S, FluxContainer, CpDgf, CmDgf> >,
  public Dune::PDELab::NumericalJacobianBoundary          <PoissonOperator<PbB, F, S, FluxContainer, CpDgf, CmDgf> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };                                // assemble boundary
  
  PoissonOperator (const PbB& pbB_, const F& f_, S& s_, FluxContainer fluxContainer_, CpDgf & cpDgf_, CmDgf & cmDgf_, unsigned int intorder_=3)  // needs boundary cond. type
    : pbB(pbB_),  f(f_), s(s_), fluxContainer(fluxContainer_), cmDgf(cmDgf_), cpDgf(cpDgf_) ,
intorder(intorder_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;
    
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;

    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV::Traits::SizeType size_type;

    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianType;
        
    // dimensions
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& 
      rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator 
           it=rule.begin(); it!=rule.end(); ++it)
      {
        
        // transform gradients from reference element to real element
        const Dune::FieldMatrix<DF,dimw,dim> 
          jac = eg.geometry().jacobianInverseTransposed(it->position());
        

        Dune::FieldVector<RF,dim> 
          globalpos = eg.geometry().global(it->position());
        RF factor = it->weight()*eg.geometry().integrationElement(it->position());
        if (s.cylindrical) 
            factor *= globalpos[1]*2*PI;

        // evaluate basis functions on reference element
        std::vector<RangeType> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        // compute u at integration point
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); i++)
          u += x(lfsu, i)*phi[i];

        Dune::FieldVector<double, 1> cp = 0;
        cpDgf.evaluate(eg.entity(), it->position(), cp);
        Dune::FieldVector<double, 1> cm = 0;
        cmDgf.evaluate(eg.entity(), it->position(), cm);
//        std::cout << globalpos << " " << cp << " " << cm << " grep" << std::endl;

        // evaluate gradient of basis functions on reference element
        std::vector<JacobianType> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);
        
        // transform gradients from reference element to real element
        std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
        for (size_type i=0; i<lfsu.size(); i++)
          jac.mv(js[i][0],gradphi[i]);

        // compute gradient of u
        Dune::FieldVector<RF,dim> gradu(0.0);
        for (size_type i=0; i<lfsu.size(); i++)
          gradu.axpy(x(lfsu, i), gradphi[i]);

        ////////////////////////// Weak Formulation /////////////////////////////////////////////
        // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
        for (size_type i=0; i<lfsu.size(); i++)
        {
          double thisResidual = ( gradu*gradphi[i] 
                                       + 1*s.l_b*4*PI*(cm-cp)*phi[i] 
                                       )*factor; 
          r.accumulate(lfsu, i, thisResidual);
        }
      }
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu, const X& x_s, 
                       const LFSV& lfsv, R& r_s) const
  {
 
    
    // some types
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;

    
    typedef typename LFSV::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;

    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV::Traits::SizeType size_type;

    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypePhi;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCp;
    typedef typename LFSU::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCm;
        

    // dimensions
    const int dim = IG::dimension;
        
    // select quadrature rule for face
    Dune::GeometryType gtface = ig.geometryInInside().type();
    const Dune::QuadratureRule<DF,dim-1>& 
      rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); 
	 it!=rule.end(); ++it)
    {
        // evaluate boundary condition type
        bool pbBtype;
        pbBtype = pbB.isDirichlet(ig,it->position());
        
        // position of quadrature point in local coordinates of element 
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
        Dune::FieldVector<DF,dim> global = ig.geometry().global(it->position());
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        if (s.cylindrical) 
            factor *= global[1]*2*PI;
        RF j;

        // evaluate basis functions on reference element
        std::vector<RangeType> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(local,phi);

        if (!pbBtype) {
 
          // evaluate flux boundary condition
          j = fluxContainer[ig.intersection().boundarySegmentIndex()][0];
              
          // integrate j
          for (size_type i=0; i<lfsv.size(); i++)
          {
            double thisResidual_s =  j*phi[i]*factor;
            r_s.accumulate(lfsu, i, thisResidual_s);       
          }
        }
        
    } 
  }

private:
  const PbB & pbB;
  const F& f;
  const S& s;
  const FluxContainer fluxContainer;
  const CpDgf & cpDgf;
  const CmDgf & cmDgf;
  unsigned int intorder;
};
