#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

/** \brief A local operator for the mass operator (L_2 integral) in the system */
template<class S>
class PnpTOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume<PnpTOperator<S> >,
    public Dune::PDELab::NumericalJacobianVolume<PnpTOperator<S> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // constructor remembers parameters
  PnpTOperator (S& s_, double tau_, unsigned int intorder_=2)
    : s(s_), tau(tau_), intorder(intorder_) {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
        
    typedef typename LFSU::template Child<0>::Type LFSU_Phi;
    typedef typename LFSU::template Child<1>::Type LFSU_Cp;
    typedef typename LFSU::template Child<2>::Type LFSU_Cm;
    const LFSU_Phi& lfsu_phi = lfsu.template getChild<0>();
    const LFSU_Cp&  lfsu_cp = lfsu.template getChild<1>();
    const LFSU_Cm&  lfsu_cm = lfsu.template getChild<2>();

    typedef typename LFSV::template Child<0>::Type LFSV_Phi;
    typedef typename LFSV::template Child<1>::Type LFSV_Cp;
    typedef typename LFSV::template Child<2>::Type LFSV_Cm;



    typedef typename LFSV_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;

    
    typedef typename LFSV_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;

    typedef typename LFSU_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV_Phi::Traits::SizeType size_type;

    // dimensions
    const int dim = EG::Geometry::dimension;

    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions on reference element
        std::vector<RangeType> phi_phi(lfsu_phi.size());
        lfsu_phi.finiteElement().localBasis().evaluateFunction(it->position(),phi_phi);
        std::vector<RangeType> phi_cp(lfsu_cp.size());
        lfsu_cp.finiteElement().localBasis().evaluateFunction(it->position(),phi_cp);
        std::vector<RangeType> phi_cm(lfsu_cm.size());
        lfsu_cm.finiteElement().localBasis().evaluateFunction(it->position(),phi_cm);

        // compute u_0, u_1 at integration point
        RF u_phi=0.0;
        for (size_type i=0; i<lfsu_phi.size(); i++)
          u_phi += x(lfsu_phi, i)*phi_phi[i];
        RF u_cp=0.0;
        for (size_type i=0; i<lfsu_cp.size(); i++)
          u_cp += x(lfsu_cp, i)*phi_cp[i];
        RF u_cm=0.0;
        for (size_type i=0; i<lfsu_cm.size(); i++)
          u_cm += x(lfsu_cm, i)*phi_cm[i];
        RF u_0=0.0;

        // integration
        Dune::FieldVector<RF,dim> 
          globalpos = eg.geometry().global(it->position());
        RF factor = it->weight() * eg.geometry().integrationElement(it->position());
        if (s.cylindrical) 
            factor *= globalpos[1]*2*PI;

        for (size_type i=0; i<lfsu_cp.size(); i++) 
          r.accumulate(lfsu_cp,i,tau*u_cp*phi_cp[i]*factor);
        for (size_type i=0; i<lfsu_cm.size(); i++) 
          r.accumulate(lfsu_cp,i,tau*u_cm*phi_cm[i]*factor);
      }
  }
private:
  double tau;
  unsigned int intorder;
  const S& s;
};
