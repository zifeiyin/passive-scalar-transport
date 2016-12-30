/*---------------------------------------------------------------------------*\
dynamicSmagorinskyST - Implementation of the dynamic Smagorinsky
		     SGS model.
    
Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010-2011 Alberto Passalacqua 
    
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynamicSmagorinskyST.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicSmagorinskyST, 0);
addToRunTimeSelectionTable(LESModel, dynamicSmagorinskyST, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dynamicSmagorinskyST::updateSubGridScaleFields
(
    const volSymmTensorField& D
)
{
    // The SGS viscosity is bounded so that nuEff cannot become negative.
    // Values are limited here, and not in nuEff, for consistency in stored
    // data and in submodels using nuSgs().
    // No warning message is printed when this limitation is applied.
    //nuSgs_ = max(cD(D),scalar(0.0))*sqr(delta())*sqrt(magSqr(D));
    nuSgs_.correctBoundaryConditions();
}

volScalarField dynamicSmagorinskyST::cD
(
    const volSymmTensorField& D
) const
{
    tmp<volSymmTensorField> LL = 
	dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    const volSymmTensorField MM
    (
        sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D))
    );

    // Locally averaging MMMM on cell faces
    volScalarField MMMM = fvc::average(magSqr(MM));

    MMMM.max(VSMALL);

    // Performing local average on cell faces on return
    return 0.5*fvc::average(LL && MM)/MMMM;
}


volScalarField dynamicSmagorinskyST::cI
(
    const volSymmTensorField& D
) const
{
    tmp<volScalarField> KK = 
	0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    const volScalarField mm
    (
        sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))))
    );

    // Locally averaging mmmm on cell faces
    volScalarField mmmm = fvc::average(magSqr(mm));

    mmmm.max(VSMALL);

    // Performing local average on cell faces on return
    return fvc::average(KK*mm)/mmmm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicSmagorinskyST::dynamicSmagorinskyST
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(typeName, U, phi, transport),
    GenEddyVisc(U, phi, transport),
    theta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "theta",
            coeffDict_,
            1.5
        )
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    invT
    (
        IOobject
        (
            "invT",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
    mesh_,
    dimensionedScalar("invT", dimVelocity/dimLength, 1)
    ),
    flm_
    (
        IOobject
        (
            "flm",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    fmm_
    (
        IOobject
        (
            "fmm",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    flm0_("flm0", flm_.dimensions(), 0.0),
    fmm0_("fmm0", fmm_.dimensions(), 1e-12),
    Prt_
    (
        IOobject
        (
            "Prt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Prt", dimLength/dimLength, 0.9)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    f1_
    (
        IOobject
        (
            "f1",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    f2_
    (
        IOobject
        (
            "f2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    T_
    (
        IOobject
        (
            "T",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    updateSubGridScaleFields(dev(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicSmagorinskyST::correct
(
    const tmp<volTensorField>& gradU
)
{
    LESModel::correct(gradU);

    const volSymmTensorField D(dev(symm(gradU)));

    k_ = cI(D)*sqr(delta())*magSqr(D);
    bound(k_,  kMin_);

    //updateSubGridScaleFields(D);

    volSymmTensorField devS(dev(symm(fvc::grad(U_))));

    volSymmTensorField Sf(dev(symm(fvc::grad(filter_(U())))));

    volSymmTensorField L(dev(filter_(sqr(U())) - (sqr(filter_(U())))));

    volSymmTensorField M(2.0*sqr(delta())*(filter_(mag(devS)*devS) - 4.0*filter_(mag(devS))*Sf));


    invT = (1.0/(theta_.value()*delta()))*pow(flm_*fmm_, 1.0/8.0);


    fvScalarMatrix flmEqn
    (
        fvm::ddt(flm_)
      + fvm::div(phi(), flm_)
     ==
        invT*(L && M)
      - fvm::Sp(invT, flm_)
    );

    flmEqn.relax();
    flmEqn.solve();

    flm_ = max(flm_, flm0_);


    fvScalarMatrix fmmEqn
    (
        fvm::ddt(fmm_)
      + fvm::div(phi(), fmm_)
     ==
        invT*(M && M)
      - fvm::Sp(invT, fmm_)
    );

    fmmEqn.relax();
    fmmEqn.solve();

    fmm_ = max(fmm_, fmm0_);

    nuSgs_ = max(flm_/fmm_,scalar(0.0))*sqr(delta())*sqrt(magSqr(D));
    nuSgs_.correctBoundaryConditions();
    volVectorField F_( filter_(U()*T_ ) - ( filter_(U())*filter_(T_) ) );
    volVectorField Tgrad_ = fvc::grad(filter_(T_));

    fvScalarMatrix f1Eqn
    (
        fvm::ddt(f1_)
      + fvm::div(phi(), f1_)
     ==
        invT*(-(L & Tgrad_) & (devS & F_ ))
      - fvm::Sp(invT, f1_)
    );

    f1Eqn.relax();
    f1Eqn.solve();

    fvScalarMatrix f2Eqn
    (
        fvm::ddt(f2_)
      + fvm::div(phi(), f2_)
     ==
        invT*(scalar(2.0) * (devS & F_ ) & (devS & F_ ))
      - fvm::Sp(invT, f2_)
    );
    f2Eqn.relax();
    f2Eqn.solve();

    Prt_ = max(f1_/f2_,scalar(0.01));

    alphat_ = nuSgs_/Prt_;
    alphat_.correctBoundaryConditions();
    

    volScalarField alphaEff("alphaEff", nu()/0.71 + alphat_);


/*
    // used for adding heat source for channel flow
    const dictionary& transportProperties = U_.db().lookupObject<IOdictionary>
    (
     "transportProperties"
    );
    dimensionedScalar targetT(transportProperties.lookup("Tbar"));

    dimensionedScalar Qh_ = (targetT-average(T_))/mag(U_.time().deltaT());
*/
    tmp<fvScalarMatrix> TEqn
    (
        fvm::ddt(T_)
      + fvm::div(phi_, T_)
      - fvm::laplacian(alphaEff, T_)
      ==
      Qh_
    );

    TEqn().relax();
    TEqn().solve();
    T_.correctBoundaryConditions();

    Info  << min(T_) << max(T_)  << average(T_) << Qh_ << endl;


}

bool dynamicSmagorinskyST::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
