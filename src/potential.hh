



template<class P, class C>
class Potential {
  public:
    P& params;
    Potential(P& _params) : params(_params) {
    }
    double operator ()(C pos) {
      const int dimworld = C::dimension;
      double dist=0;
      for (int i=0; i<dimworld; i++) {
        dist += pos[i]*pos[i];
      }
      if (dist > 1) 
        return 0;
      else
        return params.K*(1-dist*dist);
    }
};


template<class P, class C>
class Force {
  public:
    P& params;
    Force(P& _params) : params(_params) {
    }
    double operator ()(C pos) {
      const int dimworld = C::dimension;
      double dist=0;
      for (int i=0; i<dimworld; i++) {
        dist += pos[i]*pos[i];
      }
      double temp(0);
      if (dist > 1) 
        return temp*dist*params.K;
      else
        return 0*temp;
    }
};
