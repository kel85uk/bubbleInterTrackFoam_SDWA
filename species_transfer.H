Info << "Solving the species field" << endl;

DT = fluidIndicator*( interface.DTA() - interface.DTB() ) + interface.DTB(); //DT is of type volScalarField

    if
    (
        CO2.boundaryField()[interface.aPatchID()].type() 
     == fixedValueFvPatchScalarField::typeName
    )
    {
        fixedValueFvPatchScalarField& CO2l =
            refCast<fixedValueFvPatchScalarField>
            (
                CO2.boundaryField()[interface.aPatchID()]
            );
          CO2l == interface.H().value()*CO2.boundaryField()[interface.bPatchID()];
//        CO2l == 3.181e-2*CO2.boundaryField()[interface.bPatchID()];
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Concentration boundary condition "
                << "at the freeSurface patch is not "
                << fixedValueFvPatchScalarField::typeName
                << exit(FatalError);
    }
    
    if
    (
        CO2.boundaryField()[interface.bPatchID()].type() 
     == fixedGradientFvPatchScalarField::typeName
    )
    {
        fixedGradientFvPatchScalarField& CO2g =
            refCast<fixedGradientFvPatchScalarField>
            (
                CO2.boundaryField()[interface.bPatchID()]
            );

        CO2g == -interface.DTA().value()*CO2.boundaryField()[interface.aPatchID()].snGrad()/interface.DTB().value();
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Concentration boundary condition "
                << "at the freeSurfaceShadow patch is not "
                << fixedGradientFvPatchScalarField::typeName
                << exit(FatalError);
    }

fvScalarMatrix CO2Eqn
(
  fvm::ddt(CO2)
  + fvm::div(phi - fvc::meshPhi(rho,U), CO2, "div(phi,CO2)")
  //+ fvm::div(phi, CO2)
  == fvm::laplacian(DT, CO2)
);

//CO2.correctBoundaryConditions();
CO2Eqn.solve();
