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
#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpicollectivecommunication.hh>


/** \brief Writer for grid data in gnuplot format
    \tparam GridType the grid
    \tparam GridView Level- or LeafGridView
*/

template<class GridView>
class DataWriter {
    
  typedef typename GridView::Grid::ctype ctype;
  enum {dimworld = GridView::dimensionworld};
  enum {dim = GridView::dimension};

  public:

    DataWriter (const GridView & _gv) : gv(_gv), 
      communicator(gv.comm()) {}
    
    /** \brief Add cell data
        \param data An ISTL compliant vector type
        \param name associated with the data
        \param filename Filename for the output (*.dat)
    */

    template <class GFS, class DataContainer>
    void writeData(const GFS& gfs, const DataContainer& data, 
            const std::string &filename)
      {
         
        typedef double Real;
        std::ofstream out;
        if ( communicator.rank() == 0 ) {
          out.open (filename.c_str(), std::ios::out);
          out.precision(5);
          out << "This is intro" << std::endl;
          out.close();
        }
        for ( int i = 0; i < communicator.size(); i++) {
            if (i==communicator.rank()) {
              out.open (filename.c_str(), std::ios::out);
              typedef typename GridView::template Codim<0>::template Partition
                  <Dune::Interior_Partition>::Iterator LeafIterator;
          // do output
              for (LeafIterator it = 
                  gv.template begin<0,Dune::Interior_Partition>();
                  it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
              {
                  Dune::FieldVector<Real, dim> evalPos = 
                    it->geometry().center();
                  Dune::FieldVector<Real, dim> local = 
                    it->geometry().local(evalPos);
                  // construct a discrete grid function for access to solution
                  typedef Dune::PDELab::DiscreteGridFunction
                    <GFS,DataContainer> DGF;
                  const DGF udgf(gfs, data);
                  typedef typename DGF::Traits::RangeType RT;
                  RT value;
                  // evaluate the potential
                  udgf.evaluate(*it, local, value);
                  Dune::FieldVector<Real,dim> gradphi;
                  typename Dune::PDELab::DiscreteGridFunctionGradient
                    < GFS, DataContainer > grads(gfs,data);
                  grads.evaluate(*it, local, gradphi);
  

                  out << std::left << std::scientific << evalPos << "\t";
                  out << std::left << value << "\t";
                  out << std::left << gradphi << std::endl;
            }
            out.close();
          }
        communicator.barrier();
          
        }
      }
  
  private:
    const GridView &gv;

    typedef typename GridView::Traits::CollectiveCommunication CollectiveCommunication;
    const CollectiveCommunication & communicator;

};

#endif // IPBS_DATAWRITER_HH
