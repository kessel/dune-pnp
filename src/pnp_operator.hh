#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>
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
template<class PhiB, class CpB, class CmB, class F, class S, class FluxContainer>
class PnpOperator : 
  public Dune::PDELab::NumericalJacobianApplyVolume<PnpOperator<PhiB, CpB, CmB, F, S, FluxContainer> >,
  public Dune::PDELab::NumericalJacobianVolume<PnpOperator<PhiB, CpB, CmB, F, S, FluxContainer> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<PnpOperator<PhiB, CpB, CmB, F, S, FluxContainer> >,
  public Dune::PDELab::NumericalJacobianBoundary<PnpOperator<PhiB, CpB, CmB, F, S, FluxContainer> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };                                // assemble boundary
  
  PnpOperator (const PhiB& phiB_, const CpB& cpB_, const CmB& cmB_, const F& f_, S& s_, FluxContainer fluxContainer_, unsigned int intorder_=2)  // needs boundary cond. type
    : phiB(phiB_), cpB(cpB_), cmB(cmB_),  f(f_), intorder(intorder_), s(s_), fluxContainer(fluxContainer_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu_s, const X& x, const LFSV& lfsv_s, R& r) const
  {


        
 
    typedef typename LFSU::template Child<0>::Type LFSU_Phi;
    typedef typename LFSU::template Child<1>::Type LFSU_Cp;
    typedef typename LFSU::template Child<2>::Type LFSU_Cm;
    const LFSU_Phi& lfsu_phi = lfsu_s.template getChild<0>();
    const LFSU_Cp&  lfsu_cp = lfsu_s.template getChild<1>();
    const LFSU_Cm&  lfsu_cm = lfsu_s.template getChild<2>();

    typedef typename LFSV::template Child<0>::Type LFSV_Phi;
    typedef typename LFSV::template Child<1>::Type LFSV_Cp;
    typedef typename LFSV::template Child<2>::Type LFSV_Cm;
    const LFSV_Phi& lfsv_phi = lfsv_s.template getChild<0>();
    const LFSV_Cp&  lfsv_cp = lfsv_s.template getChild<1>();
    const LFSV_Cm&  lfsv_cm = lfsv_s.template getChild<2>();
    
    // some types
    typedef typename LFSV_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;

    
    typedef typename LFSV_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;

    typedef typename LFSU_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV_Phi::Traits::SizeType size_type;

    typedef typename LFSU_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypePhi;
    typedef typename LFSU_Cp::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCp;
    typedef typename LFSU_Cm::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCm;
  

        
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

        // Commons:
        
        // transform gradients from reference element to real element
        const Dune::FieldMatrix<DF,dimw,dim> 
          jac = eg.geometry().jacobianInverseTransposed(it->position());
        

        Dune::FieldVector<RF,dim> 
          globalpos = eg.geometry().global(it->position());
        RF F = 0; 
        RF a = 0; 
        RF factor = it->weight()*eg.geometry().integrationElement(it->position());

        // evaluate basis functions on reference element
        std::vector<RangeType> phi_phi(lfsu_phi.size());
        lfsu_phi.finiteElement().localBasis().evaluateFunction(it->position(),phi_phi);
        std::vector<RangeType> phi_cp(lfsu_cp.size());
        lfsu_cp.finiteElement().localBasis().evaluateFunction(it->position(),phi_cp);
        std::vector<RangeType> phi_cm(lfsu_cm.size());
        lfsu_cm.finiteElement().localBasis().evaluateFunction(it->position(),phi_cm);

        // compute u at integration point
        RF u_phi=0.0;
        for (size_type i=0; i<lfsu_phi.size(); i++)
          u_phi += x[lfsu_phi.localIndex(i)]*phi_phi[i];
        RF u_cp=0.0;
        for (size_type i=0; i<lfsu_cp.size(); i++)
          u_cp += x[lfsu_cp.localIndex(i)]*phi_cp[i];
        RF u_cm=0.0;
        for (size_type i=0; i<lfsu_cm.size(); i++)
          u_cm += x[lfsu_cm.localIndex(i)]*phi_cm[i];

        // evaluate gradient of basis functions on reference element
        std::vector<JacobianTypePhi> js_phi(lfsu_phi.size());
        lfsu_phi.finiteElement().localBasis().evaluateJacobian(it->position(),js_phi);
        std::vector<JacobianTypeCp> js_cp(lfsu_cp.size());
        lfsu_cp.finiteElement().localBasis().evaluateJacobian(it->position(),js_cp);
        std::vector<JacobianTypeCp> js_cm(lfsu_cm.size());
        lfsu_cm.finiteElement().localBasis().evaluateJacobian(it->position(),js_cm);
        
        // transform gradients from reference element to real element
        std::vector<Dune::FieldVector<RF,dim> > gradphi_phi(lfsu_phi.size());
        for (size_type i=0; i<lfsu_phi.size(); i++)
          jac.mv(js_phi[i][0],gradphi_phi[i]);
        std::vector<Dune::FieldVector<RF,dim> > gradphi_cp(lfsu_cp.size());
        for (size_type i=0; i<lfsu_cp.size(); i++)
          jac.mv(js_cp[i][0],gradphi_cp[i]);
        std::vector<Dune::FieldVector<RF,dim> > gradphi_cm(lfsu_cm.size());
        for (size_type i=0; i<lfsu_cm.size(); i++)
          jac.mv(js_cm[i][0],gradphi_cm[i]);

        // compute gradient of u
        Dune::FieldVector<RF,dim> gradu_phi(0.0);
        for (size_type i=0; i<lfsu_phi.size(); i++)
          gradu_phi.axpy(x[lfsu_phi.localIndex(i)],gradphi_phi[i]);
        Dune::FieldVector<RF,dim> gradu_cp(0.0);
        for (size_type i=0; i<lfsu_cp.size(); i++)
          gradu_cp.axpy(x[lfsu_cp.localIndex(i)],gradphi_cp[i]);
        Dune::FieldVector<RF,dim> gradu_cm(0.0);
        for (size_type i=0; i<lfsu_cm.size(); i++)
          gradu_cm.axpy(x[lfsu_cm.localIndex(i)],gradphi_cm[i]);

        ////////////////////////// Weak Formulation /////////////////////////////////////////////
        // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
        for (size_type i=0; i<lfsu_phi.size(); i++)
          r[lfsu_phi.localIndex(i)] += ( gradu_phi*gradphi_phi[i] + a*u_phi*phi_phi[i] - F*phi_phi[i] )*factor; 

        ////////////////////////// NOW CP /////////////////////////////////////////////
        
        // integrate grad u * grad cp_i + a*u*cp_i - f cp_i
        for (size_type i=0; i<lfsu_cp.size(); i++)
          r[lfsu_cp.localIndex(i)] += ( gradu_cp*gradphi_cp[i] + a*u_cp*phi_cp[i] - F*phi_cp[i] )*factor; 

        ///////////////////////// Finally CM ////////////////////////////////////////

        // integrate grad u * grad cm_i + a*u*cm_i - f cm_i
        for (size_type i=0; i<lfsu_cm.size(); i++)
          r[lfsu_cm.localIndex(i)] += ( gradu_cm*gradphi_cm[i] + a*u_cm*phi_cm[i] - F*phi_cm[i] )*factor; 
      }
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s, 
                       const LFSV& lfsv_s, R& r_s) const
  {
 
    typedef typename LFSU::template Child<0>::Type LFSU_Phi;
    typedef typename LFSU::template Child<1>::Type LFSU_Cp;
    typedef typename LFSU::template Child<2>::Type LFSU_Cm;
    const LFSU_Phi& lfsu_phi = lfsu_s.template getChild<0>();
    const LFSU_Cp&  lfsu_cp = lfsu_s.template getChild<1>();
    const LFSU_Cm&  lfsu_cm = lfsu_s.template getChild<2>();

    typedef typename LFSV::template Child<0>::Type LFSV_Phi;
    typedef typename LFSV::template Child<1>::Type LFSV_Cp;
    typedef typename LFSV::template Child<2>::Type LFSV_Cm;
    const LFSV_Phi& lfsv_phi = lfsv_s.template getChild<0>();
    const LFSV_Cp&  lfsv_cp = lfsv_s.template getChild<1>();
    const LFSV_Cm&  lfsv_cm = lfsv_s.template getChild<2>();
    
    // some types
    typedef typename LFSV_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::DomainFieldType DF;

    
    typedef typename LFSV_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeFieldType RF;

    typedef typename LFSU_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::RangeType RangeType;

    typedef typename LFSV_Phi::Traits::SizeType size_type;

    typedef typename LFSU_Phi::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypePhi;
    typedef typename LFSU_Cp::Traits::FiniteElementType::
      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCp;
    typedef typename LFSU_Cm::Traits::FiniteElementType::
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
        typename PhiB::Traits::RangeType phiBtype;
        phiB.evaluate(ig,it->position(),phiBtype);
        
        typename CpB::Traits::RangeType cpBtype;
        cpB.evaluate(ig,it->position(),cpBtype);
        
        typename CpB::Traits::RangeType cmBtype;
        cmB.evaluate(ig,it->position(),cmBtype);
          
        // position of quadrature point in local coordinates of element 
        Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());
        RF j;


        // evaluate basis functions on reference element
        std::vector<RangeType> phi_phi(lfsu_phi.size());
        lfsu_phi.finiteElement().localBasis().evaluateFunction(local,phi_phi);
        std::vector<RangeType> phi_cp(lfsu_cp.size());
        lfsu_cp.finiteElement().localBasis().evaluateFunction(local,phi_cp);
        std::vector<RangeType> phi_cm(lfsu_cm.size());
        lfsu_cm.finiteElement().localBasis().evaluateFunction(local,phi_cm);

        if (phiBtype < 1) {
 
          // evaluate flux boundary condition
          j = fluxContainer[ig.intersection().boundarySegmentIndex()][0];
              
          // integrate j
          for (size_type i=0; i<lfsv_phi.size(); i++)
            r_s[lfsu_phi.localIndex(i)] += j*phi_phi[i]*factor;
        }
        
        if (cpBtype < 1) {
          // evaluate flux boundary condition
          j = fluxContainer[ig.intersection().boundarySegmentIndex()][1];
 
          // integrate j
          for (size_type i=0; i<lfsv_cp.size(); i++)
            r_s[lfsu_cp.localIndex(i)] += j*phi_cp[i]*factor;
        }
        if (cmBtype < 1) {
              
          // evaluate flux boundary condition
          j = fluxContainer[ig.intersection().boundarySegmentIndex()][2];
              
          // integrate j
          for (size_type i=0; i<lfsv_cm.size(); i++)
            r_s[lfsu_cm.localIndex(i)] += j*phi_cm[i]*factor;
        }
    }

  }

private:
  const PhiB & phiB;
  const CpB & cpB;
  const CmB & cmB;
  const F& f;
  const S& s;
  const FluxContainer fluxContainer;
  unsigned int intorder;
};
