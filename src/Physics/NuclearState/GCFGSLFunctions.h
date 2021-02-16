

namespace genie{
  namespace gcf{
    using spline_ptr = std::shared_ptr<Spline>;
    class PDF_Pr: public ROOT::Math::IBaseFunctionOneDim{
      public: 
        PDF_Pr(const spline_ptr &p):IBaseFunctionOneDim(),fSpline(p){;}
        PDF_Pr() = delete;
        
        // ROOT::Math::IBaseFunctionOneDim interface
        unsigned int                      NDim   (void)             const{return 1;}
        double                            DoEval (double xin) const{return fSpline->Evaluate(xin);}
        ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const{return new PDF_Pr(*this);}


      private:
        spline_ptr fSpline;
    };

  }
}
