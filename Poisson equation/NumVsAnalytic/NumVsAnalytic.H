#ifndef NumVsAnalytic_H
#define NumVsAnalytic_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class NumVsAnalytic Declaration
\*---------------------------------------------------------------------------*/

class NumVsAnalytic
:
    public fvMeshFunctionObject
{

	const Time& time_; 
        const fvMesh& mesh_; 

        volScalarField phi_;
        volScalarField phiAnalytic_;
        volScalarField err_;

public:

    //- Runtime type information
    TypeName("NumVsAnalytic");


    // Constructors

        //- Construct from Time and dictionary
        NumVsAnalytic
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        NumVsAnalytic(const NumVsAnalytic&) = delete;

        //- No copy assignment
        void operator=(const NumVsAnalytic&) = delete;

	const Time& getTime() const
        {
            return time_; 
        }

    //- Destructor
    virtual ~NumVsAnalytic() = default;

   
   //- Functions
   virtual bool execute();   //- Execute, currently does nothing
   virtual bool start();
   virtual bool end();   //- Execute at the final time-loop, currently does nothing  
   virtual bool write();     //- Write the NumVsAnalytic
   virtual void report();
   
   void updateMesh(const mapPolyMesh&);
   void movePoints(const polyMesh&);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
