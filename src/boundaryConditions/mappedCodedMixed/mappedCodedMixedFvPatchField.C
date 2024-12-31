/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*--------------------------------------------------------------------------*/

#include "mappedCodedMixedFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedCodedMixedFvPatchField<Type>::mappedCodedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
  //fixedValueFvPatchField<Type>(p, iF),
    codedMixedFvPatchField<Type>(p, iF),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this),
    window_(1.0)
{}


template<class Type>
Foam::mappedCodedMixedFvPatchField<Type>::mappedCodedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
  //fixedValueFvPatchField<Type>(p, iF, dict),
    codedMixedFvPatchField<Type>(p, iF, dict),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, dict),
    window_(readScalar(dict.lookup("window")))
{}


template<class Type>
Foam::mappedCodedMixedFvPatchField<Type>::mappedCodedMixedFvPatchField
(
    const mappedCodedMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
  //fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    codedMixedFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, ptf),
    window_(ptf.window_)
{}


template<class Type>
Foam::mappedCodedMixedFvPatchField<Type>::mappedCodedMixedFvPatchField
(
    const mappedCodedMixedFvPatchField<Type>& ptf
)
:
  //fixedValueFvPatchField<Type>(ptf),
    codedMixedFvPatchField<Type>(ptf),
    mappedPatchFieldBase<Type>(ptf),
    window_(ptf.window_)
{}


template<class Type>
Foam::mappedCodedMixedFvPatchField<Type>::mappedCodedMixedFvPatchField
(
    const mappedCodedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
  //fixedValueFvPatchField<Type>(ptf, iF),
    codedMixedFvPatchField<Type>(ptf, iF),
    mappedPatchFieldBase<Type>(this->mapper(this->patch(), iF), *this, ptf),
    window_(ptf.window_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::mappedPatchBase& Foam::mappedCodedMixedFvPatchField<Type>::mapper
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


template<class Type>
void Foam::mappedCodedMixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

   //This part set the inflow value to the sampled value,
  //but we don't want that.  Now this BC should act just
  //like a codedMixed BC until we modify it, but we at
  //least have access to the sampled field now.
  //this->operator==(this->mappedField());

    if (debug)
    {
        Info<< "mapped on field:"
            << this->internalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

  //fixedValueFvPatchField<Type>::updateCoeffs();
    codedMixedFvPatchField<Type>::updateCoeffs(); // is this where the code part of the BC get evaluated?

    // Grab the upstream information
    Field<Type> sampledField = this->mappedField();

    //Initialize the sampledFieldAvg with the first sampled value
    if (this->db().time().value() == this->db().time().startTime().value())
    {
        sampledFieldAvg = sampledField;    
    }
    else
    {
        //Begin calculating the backward time-average
        scalar dt = this->db().time().deltaTValue();
        scalar beta = dt/window_;
        sampledFieldAvg = (1 - beta)*sampledFieldAvg + beta*sampledField;
    }

    Field<Type> perturbation = sampledField - sampledFieldAvg;

    this->refValue() += perturbation;
    
}


template<class Type>
void Foam::mappedCodedMixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchFieldBase<Type>::write(os);
    codedMixedFvPatchField<Type>::write(os);
    os.writeKeyword("window") << window_ << token::END_STATEMENT << nl;
 
  //this->refValue().writeEntry("refValue", os);
}


// ************************************************************************* //
