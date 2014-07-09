/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM 
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Fabian K.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Wale.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Wale, 0);
addToRunTimeSelectionTable(LESModel, Wale, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Wale::updateSubGridScaleFields(const volTensorField& gradU)
{
    nuSgs_ = ck_*delta()*sqrt(k(gradU));
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Wale::Wale
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

    cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cw",
            coeffDict_,
            0.1
        )
     ),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),

    smallConst_("smallConst_",dimensionSet(0, 0, -10, 0, 0),SMALL)

{
    updateSubGridScaleFields(fvc::grad(U));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Wale::correct(const tmp<volTensorField>& gradU)
{
    GenEddyVisc::correct(gradU);
    updateSubGridScaleFields(gradU());
}


bool Wale::read()
{
    if (GenEddyVisc::read())
    {
        cw_.readIfPresent(coeffDict());

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
