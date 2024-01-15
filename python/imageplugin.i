%module imageplugin

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
%include "std_pair.i"
%include "std_string.i"

namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "OpenMMImage.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"

using namespace OpenMM;

%}

%pythoncode %{
import openmm as mm
import openmm.unit as unit
%}

/*
 * Add units to function outputs.
*/
%pythonappend OpenMM::ImageLangevinIntegrator::getTemperature() const %{
    val=unit.Quantity(val, unit.kelvin)
%}

%pythonappend OpenMM::ImageLangevinIntegrator::getFriction() const %{
    val=unit.Quantity(val, 1/unit.picosecond)
%}

// %pythonappend OpenMM::ImageLangevinIntegrator::getCellSize() const %{
//     val=unit.Quantity(val, unit.nanometer)
// %}


%pythonappend OpenMM::ImageIntegrator::getCellSize() const %{
    val=unit.Quantity(val, unit.nanometer)
%}

// 
%pythonappend OpenMM::ImageCustomIntegrator::addTabulatedFunction(const std::string & name, TabulatedFunction * function) %{
   function.thisown=0
%}

%pythonprepend OpenMM::ImageCustomIntegrator::addTabulatedFunction(const std::string &name, TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyObject* image = PyImport_AddModule("imageplugin");
        PyObject* image_exception = PyObject_GetAttrString(image, "OpenMMException");
        PyErr_SetString(image_exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

%inline %{
    #include <cstring>
    #include <numpy/arrayobject.h>
%}
%header %{
namespace OpenMM {

PyObject *copyVVec3ToList(std::vector<Vec3> vVec3) {
  int i, n;
  PyObject *pyList;

  n=vVec3.size(); 
  pyList=PyList_New(n);
  // get the attribute Vec3 in openmm module
  PyObject* mm = PyImport_AddModule("openmm");
  PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
  for (i=0; i<n; i++) {
    OpenMM::Vec3& v = vVec3.at(i);
    PyObject* args = Py_BuildValue("(d,d,d)", v[0], v[1], v[2]);
    PyObject* pyVec = PyObject_CallObject(vec3, args);
    Py_DECREF(args);
    PyList_SET_ITEM(pyList, i, pyVec);
  }
  return pyList;
}

int isNumpyAvailable() {
    static bool initialized = false;
    static bool available = false;
    if (!initialized) {
        initialized = true;
        available = (_import_array() >= 0);
    }
    return available;
}

} // namespace OpenMM
%}

%extend OpenMM::ImageCustomIntegrator {
    PyObject* getPerDofVariable(int index) const {
        std::vector<Vec3> values;
        self->getPerDofVariable(index, values);
        return copyVVec3ToList(values);
    }
}

namespace OpenMM {

class ImageIntegrator {
public:
    ImageIntegrator(int numCells=2, double zmax=-1);

    int getNumCells() const;
    void setNumCells(int cells);
    double getCellSize() const;
    void setCellSize(double zmax);
    void setUseImageParticle(bool use);
    const bool& getUseImageParticle() const;
    int setImagePair(int image, int parent);
    const std::vector<std::pair<int,int> >& getImagePairs() const;
    bool isParticleImage(int i) const;
};

class ImageLangevinIntegrator : public Integrator, public ImageIntegrator {
public:
    ImageLangevinIntegrator(double temperature, double frictionCoeff, double stepSize, int numCells=2, double zmax=-1);

    double getTemperature() const;
    void setTemperature(double temp);
    double getFriction() const;
    void setFriction(double coeff);
    int getRandomNumberSeed() const;
    void setRandomNumberSeed(int seed);
    virtual void step(int steps);
};

class ImageCustomIntegrator : public Integrator, public ImageIntegrator {
public:
   enum ComputationType {
      ComputeGlobal =  0,
      ComputePerDof =  1,
      ComputeSum =  2,
      ConstrainPositions =  3,
      ConstrainVelocities =  4,
      UpdateContextState =  5,
      IfBlockStart =  6,
      WhileBlockStart =  7,
      BlockEnd =  8
   };

   ImageCustomIntegrator(double stepSize, int numCells=2, double zmax=-1);
   ~ImageCustomIntegrator();

   int getNumGlobalVariables() const ;
   int getNumPerDofVariables() const ;
   int getNumComputations() const ;
   int getNumTabulatedFunctions() const ;
   int addGlobalVariable(const std::string &name, double initialValue) ;
   const std::string& getGlobalVariableName(int index) const ;
   int addPerDofVariable(const std::string &name, double initialValue) ;
   const std::string& getPerDofVariableName(int index) const ;
   double getGlobalVariable(int index) const ;
   double getGlobalVariableByName(const std::string &name) const ;
   void setGlobalVariable(int index, double value) ;
   void setGlobalVariableByName(const std::string &name, double value) ;
   %apply std::vector< Vec3 > & OUTPUT { std::vector< Vec3 > & values };
   void getPerDofVariable(int index, std::vector< Vec3 > &values) const ;
   %clear std::vector< Vec3 > & values;
   %apply std::vector< Vec3 > & OUTPUT { std::vector< Vec3 > & values };
   void getPerDofVariableByName(const std::string &name, std::vector< Vec3 > &values) const ;
   %clear std::vector< Vec3 > & values;
   void setPerDofVariable(int index, const std::vector< Vec3 > &values) ;
   void setPerDofVariableByName(const std::string &name, const std::vector< Vec3 > &values) ;
   int addComputeGlobal(const std::string &variable, const std::string &expression) ;
   int addComputePerDof(const std::string &variable, const std::string &expression) ;
   int addComputeSum(const std::string &variable, const std::string &expression) ;
   int addConstrainPositions() ;
   int addConstrainVelocities() ;
   int addUpdateContextState() ;
   int beginIfBlock(const std::string &condition) ;
   int beginWhileBlock(const std::string &condition) ;
   int endBlock() ;
   %apply int & OUTPUT { ComputationType & type };
   %apply std::string & OUTPUT { std::string & variable };
   %apply std::string & OUTPUT { std::string & expression };
   void getComputationStep(int index, ComputationType &type, std::string &variable, std::string &expression) const ;
   %clear ComputationType & type;
   %clear std::string & variable;
   %clear std::string & expression;
   int addTabulatedFunction(const std::string &name, TabulatedFunction *function) ;
   const TabulatedFunction& getTabulatedFunction(int index) const ;
   TabulatedFunction& getTabulatedFunction(int index) ;
   const std::string& getTabulatedFunctionName(int index) const ;
   const std::string& getKineticEnergyExpression() const ;
   void setKineticEnergyExpression(const std::string &expression) ;
   int getRandomNumberSeed() const ;
   void setRandomNumberSeed(int seed) ;
   virtual void step(int steps) ;
};


}
