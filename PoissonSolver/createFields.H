Info<< "Reading field phi\n" << endl;
volScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar 
    (
        "phi",     
        dimensionSet(0, 2, 0, 0, 0, 0, 0), 
        0.0       
    )
);

volScalarField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
