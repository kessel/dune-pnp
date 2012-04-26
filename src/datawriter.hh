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

    template <class GFS, class DataContainer>
    void writeIpbsCellData(const GFS& gfs, const DataContainer& data, 
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
    
        // do output
        for (LeafIterator it = 
            gv.template begin<codim,Dune::Interior_Partition>();
            it!=gv.template end<codim,Dune::Interior_Partition>(); ++it)
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
            out << std::left << gradphi << "\t";
            out << std::left << gradphi.two_norm() 
                  << "\n";
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
