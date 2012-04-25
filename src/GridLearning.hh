
#include <iostream>


template<class GV>
class GridLearning {
  public:
    GridLearning(GV& gv);
    ~GridLearning();
    void coutElementPositions();
    void coutCornerPositions();
    void coutCornerByElement();
  private:
    typedef GV myGV;
    myGV& gv;
};



//template<class GV>
//GridLearning::GridLearning(GV &gv) {
//  self->gv = gv;
//}



template<class GV>
GridLearning<GV>::GridLearning(GV & _gv) : gv(_gv) {
}


template<class GV>
GridLearning<GV>::~GridLearning() {
}

template<class GV>
void GridLearning<GV>::coutElementPositions() {

  typedef typename myGV::template Codim<0>::Iterator GVIT;
  
  GVIT it = gv.template begin<0>();

  for (; it != gv.template end<0>(); ++it) {
    std::cout << it->geometry().center() << std::endl;
  }
}

template<class GV>
void GridLearning<GV>::coutCornerPositions() {
  const int dimworld = GV::dimensionworld;
  
  typedef typename myGV::template Codim<dimworld>::Iterator CornerIterator;

  CornerIterator it = gv.template begin<dimworld>() ;

  for (; it != gv.template end<dimworld>(); ++it) {
    std::cout << "corner " << it->geometry().center() << std::endl;
  }
  
}



template<class GV>
void GridLearning<GV>::coutCornerByElement() {

  typedef typename myGV::template Codim<0>::Iterator GVIT;
  typedef typename myGV::IntersectionIterator IntersectionIterator;
  
  GVIT it = gv.template begin<0>();

  for (; it != gv.template end<0>(); ++it) {
    std::cout << "intersections: ";
    for (IntersectionIterator iit = gv.ibegin(*it); iit!=  gv.iend(*it); ++iit) {
      std::cout << " (" << iit->geometry().center() << ") " ;
    }
    std::cout << std::endl;

  }
}
