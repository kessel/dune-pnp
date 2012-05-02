#ifndef IPBS_DATAWRITER_HH
#define IPBS_DATAWRITER_HH

/** @file
    @author Alexander Schlaich
    @brief Provides gnuplot readable output
    based on the Dune::Gnuplotwriter,
    but iterating over elements for element-data
*/

#include <string>
#include <iostream>
#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

/** \brief Writer for grid data in gnuplot format
    \tparam GridType the grid
    \tparam GridView Level- or LeafGridView
*/

template<class GridView, class Real, int codim=0>
class DataWriter {
    
  typedef typename GridView::Grid::ctype ctype;
  enum {dimworld = GridView::dimensionworld};
  enum {dim = GridView::dimension};

  public:
    DataWriter (const GridView & _gv, Dune::MPIHelper & _helper)
      : gv(_gv), helper(_helper)
      {}
    
    /** \brief Add cell data
        \param data An ISTL compliant vector type
        \param name associated with the data
        \param filename Filename for the output (*.dat)
    */

    template <class GFS, class U>
    void writePNPCellData(GFS& gfs, const U& u, 
            const std::string& name, const std::string &filename)
      {
        std::string myFilename = filename_helper(filename);
        std::ofstream out;
        out.open (myFilename.c_str(), std::ios::out);
        out.precision(5);

        out << std::left << std::setw(dim*12+1) << "#coordinates\t";
        out << std::left << std::setw(12) << "potential\t";
        out << std::left << std::setw(dim*12+1) << "electric field\t";
        out << std::left << std::setw(12) << "|E|\t" << "total ion density";
        out << std::endl << std::endl;

        typedef typename GridView::template Codim<codim>::template Partition
            <Dune::Interior_Partition>::Iterator LeafIterator;

//    typedef typename LFSU::template Child<0>::Type LFSU_Phi;
//    typedef typename LFSU::template Child<1>::Type LFSU_Cp;
//    typedef typename LFSU::template Child<2>::Type LFSU_Cm;
//    LFSU_Phi& lfsu_phi(lfsu_s.template getChild<0>());
//    LFSU_Cp&  lfsu_cp(lfsu_s.template getChild<1>());
//    LFSU_Cm&  lfsu_cm(lfsu_s.template getChild<2>());
//
//
//    // some types
//    typedef typename LFSU_Phi::Traits::FiniteElementType::
//      Traits::LocalBasisType::Traits::DomainFieldType DF;
//
//    typedef typename LFSU_Phi::Traits::FiniteElementType::
//      Traits::LocalBasisType::Traits::RangeFieldType RF;
//
//    typedef typename LFSU_Phi::Traits::FiniteElementType::
//      Traits::LocalBasisType::Traits::RangeType RangeType;
//
//    typedef typename LFSU_Phi::Traits::SizeType size_type;
//
//    typedef typename LFSU_Phi::Traits::FiniteElementType::
//      Traits::LocalBasisType::Traits::JacobianType JacobianTypePhi;
//    typedef typename LFSU_Cp::Traits::FiniteElementType::
//      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCp;
//    typedef typename LFSU_Cm::Traits::FiniteElementType::
//      Traits::LocalBasisType::Traits::JacobianType JacobianTypeCm;
//  
//        
//    // dimensions
//    const int dim = GridView::dimension;
//    const int dimw = GridView::dimensionworld;
//
//    // select quadrature rule
//
//    // loop over quadrature points
        for (LeafIterator it = 
            gv.template begin<codim,Dune::Interior_Partition>();
            it!=gv.template end<codim,Dune::Interior_Partition>(); ++it)
      {
//    lfsu_phi.bind(*it);
//    lfsu_cp.bind(*it);
//    lfsu_cm.bind(*it);
//           Dune::FieldVector<Real, dim> local = 
//              it->geometry().local(it->geometry().center());
//
//        // Commons:
//        
//        // transform gradients from reference element to real element
//        const Dune::FieldMatrix<DF,dimw,dim> 
//          jac = it->geometry().jacobianInverseTransposed(local);
//        
//
//        Dune::FieldVector<RF,dim> 
//          globalpos = it->geometry().global(local);
//        RF F = 0; 
//        RF a = 0; 
////        RF factor = it->weight()*it->geometry().integrationElement(local);
//
//        // evaluate basis functions on reference element
//        std::vector<RangeType> phi_phi(lfsu_phi.size());
//        lfsu_phi.finiteElement().localBasis().evaluateFunction(local,phi_phi);
//        std::vector<RangeType> phi_cp(lfsu_cp.size());
//        lfsu_cp.finiteElement().localBasis().evaluateFunction(local,phi_cp);
//        std::vector<RangeType> phi_cm(lfsu_cm.size());
//        lfsu_cm.finiteElement().localBasis().evaluateFunction(local,phi_cm);
//
//        // compute u at integration point
//        RF u_phi=0.0;
//        for (size_type i=0; i<lfsu_phi.size(); i++)
//          u_phi += x[lfsu_phi.localIndex(i)]*phi_phi[i];
//        RF u_cp=0.0;
//        for (size_type i=0; i<lfsu_cp.size(); i++)
//          u_cp += x[lfsu_cp.localIndex(i)]*phi_cp[i];
//        RF u_cm=0.0;
//        for (size_type i=0; i<lfsu_cm.size(); i++)
//          u_cm += x[lfsu_cm.localIndex(i)]*phi_cm[i];

        // transform gradients from reference element to real element
//        std::vector<Dune::FieldVector<RF,dim> > gradphi_phi(lfsu_phi.size());
//        for (size_type i=0; i<lfsu_phi.size(); i++)
//          jac.mv(js_phi[i][0],gradphi_phi[i]);
//        std::vector<Dune::FieldVector<RF,dim> > gradphi_cp(lfsu_cp.size());
//        for (size_type i=0; i<lfsu_cp.size(); i++)
//          jac.mv(js_cp[i][0],gradphi_cp[i]);
//        std::vector<Dune::FieldVector<RF,dim> > gradphi_cm(lfsu_cm.size());
//        for (size_type i=0; i<lfsu_cm.size(); i++)
//          jac.mv(js_cm[i][0],gradphi_cm[i]);
//
//        // compute gradient of u
//        Dune::FieldVector<RF,dim> gradu_phi(0.0);
//        for (size_type i=0; i<lfsu_phi.size(); i++)
//          gradu_phi.axpy(x[lfsu_phi.localIndex(i)],gradphi_phi[i]);
//        Dune::FieldVector<RF,dim> gradu_cp(0.0);
//        for (size_type i=0; i<lfsu_cp.size(); i++)
//          gradu_cp.axpy(x[lfsu_cp.localIndex(i)],gradphi_cp[i]);
//        Dune::FieldVector<RF,dim> gradu_cm(0.0);
//        for (size_type i=0; i<lfsu_cm.size(); i++)
//          gradu_cm.axpy(x[lfsu_cm.localIndex(i)],gradphi_cm[i]);

 
    
//        // do output
            Dune::FieldVector<Real, dim> evalPos = 
              it->geometry().center();
            Dune::FieldVector<Real, dim> local = 
              it->geometry().local(evalPos);

  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,0> PhiGFSS;
  PhiGFSS phiGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<PhiGFSS, U> PhiDGF;
  PhiDGF phiDGF(phiGFSS, u);
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,1> CpGFSS;
  CpGFSS cpGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<CpGFSS, U> CpDGF;
  CpDGF cpDGF(cpGFSS, u);
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS,2> CmGFSS;
  CmGFSS cmGFSS(gfs);
  typedef Dune::PDELab::DiscreteGridFunction<CmGFSS, U> CmDGF;
  CmDGF cmDGF(cmGFSS, u);

            typedef typename PhiDGF::Traits::RangeType phiRT;
            phiRT phi;
            phiDGF.evaluate(*it, local, phi);
            typedef typename PhiDGF::Traits::RangeType cpRT;
            phiRT cp;
            cpDGF.evaluate(*it, local, cp);
            typedef typename PhiDGF::Traits::RangeType cmRT;
            phiRT cm;
            cmDGF.evaluate(*it, local, cm);


//            Dune::FieldVector<Real,dim> gradphi;
//            typename Dune::PDELab::DiscreteGridFunctionGradient
//              < GFS, DataContainer > grads(gfs,data);
//            grads.evaluate(*it, local, gradphi);

            out << std::left << std::scientific << it->geometry().center() << " "
              << phi << " " << cp << " " << cm << std::endl;
//            out << std::left << u_phi << "\t";
//            out << std::left << u_cp << "\t";
//            out << std::left << u_cm 
//                  << "\n";
        }
        
        out.close();
      }
  
  private:
    const GridView &gv;
    const Dune::MPIHelper &helper;

    std::string filename_helper(const std::string &name)
    {
      // generate filename for process data
      std::ostringstream pieceName;
      if( helper.size() > 1 )
      {
        pieceName << "s" << std::setfill( '0' ) << std::setw( 4 ) << helper.size() << ":";
        pieceName << "p" << std::setfill( '0' ) << std::setw( 4 ) << helper.rank() << ":";
      }
      pieceName << name << ".dat";
      return pieceName.str();
    }

};

#endif // IPBS_DATAWRITER_HH
