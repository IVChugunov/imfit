#include "function_object.h"



/// \brief Class for image function using analytic edge-on exponential disk and
///        vertical generalized-secant function
class EdgeOnBrokenDisk : public FunctionObject
{
  // the following static constant will be defined/initialized in the .cpp file
  static const char  className[];
  
  public:
    // Constructors:
    EdgeOnBrokenDisk( );
    // redefined method/member function:
    void  Setup( double params[], int offsetIndex, double xc, double yc );
    double  GetValue( double x, double y );
    // No destructor for now

    // class method for returning official short name of class
    static void GetClassShortName( string& classname ) { classname = className; };


  protected:
    double CalculateIntensity( double r, double z );
    int  CalculateSubsamples( double r, double z );


  private:
    double IntegrandExp(double u, void *p);
    double Luminosity(double R);
    double LOSIntegral(double r);
    double RadialValue(double r);
    void   BuildRadialTable();
    double  x0, y0, PA, L_0, h1, h2, r_break, n, z_0;   // parameters
    int TABLE_SIZE;
    double  PA_rad, cosPA, sinPA;   // other useful quantities
    vector<double> radialTable;
    double  alpha, scaledZ0, Sigma_00_2, two_to_alpha;   // other useful quantities
};
