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

Class
    Foam::incompressible::LESModels::Wale

Description
    The Wale Model for incompressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        D = symm(grad(U));        
        k = nuSgs^2/(ck^2*delta^2);

   ->   nuSgs calculation and therefor k calculation like
        LES-Simulation for Acoustics (Wagner/Huttl/Sagaut), page 116
   ->   cw in this implementation = cw^2

   -> with the definition of nuSgs k is calcutated as 
      (: correspond to the double inner product and * to the inner product)

        k = (cw*delta/ck)^2*(Sd:Sd)^3/((D:D)^(5/2)+(Sd:Sd)^(5/4))^2

   -> calculation of nuSgs and nuEff
 
        nuSgs = ck*sqrt(k)*delta
        nuEff = nuSgs + nu

   -> calculation of Sd (: correspond to the double inner product 
                           and * to the inner product)

       Sd = symm(grad(U))*symm(grad(U)) + skew(grad(U))*skew(grad(U)) - 
          - 1./3.*I*(symm(grad(U)):symm(grad(U)) - skew(grad(U)):skew(grad(U)))

    \endverbatim

SourceFiles
    Wale.C

\*---------------------------------------------------------------------------*/

#ifndef Wale_H
#define Wale_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class Wale Declaration
\*---------------------------------------------------------------------------*/

class Wale : public GenEddyVisc
{
    // Private data

        dimensionedScalar cw_;
        dimensionedScalar ck_;

        dimensionedScalar smallConst_;
    // 
        tmp<volTensorField> Sd(const tmp<volTensorField>& gradU) const
        {
	  return ((symm(gradU)&symm(gradU))+(skew(gradU)&skew(gradU))-
		  ((1.0/3.0)*I)*((symm(gradU)&&symm(gradU))-
                  (skew(gradU)&&skew(gradU))));
        }

    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volTensorField& gradU);

        // Disallow default bitwise copy construct and assignment
        Wale(const Wale&);
        Wale& operator=(const Wale&);


public:

    //- Runtime type information
    TypeName("Wale");

    // Constructors

        //- Construct from components
        Wale
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~Wale()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const
        {
          return 
	  (	 
	      sqr(cw_*delta()/ck_)*
	      (  
	         pow(Sd(gradU)&&Sd(gradU),3.0)/
	         (
		     pow(pow((symm(gradU)&&symm(gradU)),5.0/2.0) + 
		     pow((Sd(gradU)&&Sd(gradU)),5.0/4.0),2.0)+smallConst_
		 )
              )
	  );
        }

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(U()));
        }


        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
