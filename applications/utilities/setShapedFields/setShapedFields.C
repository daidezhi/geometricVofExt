/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Dezhi Dai, Argonne National Laboratory (ANL)
-------------------------------------------------------------------------------
License
    This file is part of geometricVofExt, which is a geometric VOF extension
    to OpenFOAM.

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

Application
    setShapedFields

Description
    Set field values using arbitrary shapes through a dictionary. For cells
    that are partially enclosed by these shapes, the values are averaged
    based on volume fractions.

\*---------------------------------------------------------------------------*/

#include "CGALSurfaceMesh.H"
#include "surfaceMeshBooleanOperation.H"
#include "cellToSurfaceMeshLocation.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "dynamicOversetFvMesh.H"
//#include "reconstruction.H"
#include "OSspecific.H"
#include <omp.h>
#include <variant>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Simple tuple of field type and name (read from stream)
class fieldDescription
{
    word type_;
    word name_;

public:

    const word& type() const noexcept { return type_; }
    const word& name() const noexcept { return name_; }

    explicit fieldDescription(Istream& is)
    {
        is >> type_;
        is >> name_;

        // Eg, read as "volScalarFieldValue", but change to "volScalarField"
        if (type_.ends_with("Value"))
        {
            type_.erase(type_.size()-5);
        }
    }
};


// Consume unused field information
template<class Type>
bool consumeUnusedType(const fieldDescription& fieldDesc, Istream& is)
{
    typedef GeometricField<Type, fvPatchField, volMesh>  fieldType;

    if (fieldDesc.type() == fieldType::typeName)
    {
        (void) pTraits<Type>(is);
        return true;
    }

    return false;
}


// Consume unused field information
static bool consumeUnused(const fieldDescription& fieldDesc, Istream& is)
{
    return
    (
        consumeUnusedType<scalar>(fieldDesc, is)
     || consumeUnusedType<vector>(fieldDesc, is)
     || consumeUnusedType<sphericalTensor>(fieldDesc, is)
     || consumeUnusedType<symmTensor>(fieldDesc, is)
     || consumeUnusedType<tensor>(fieldDesc, is)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Setting volume fields
template<class Type>
bool setCellFieldType
(
    const fieldDescription& fieldDesc,
    const fvMesh& mesh,
    const DynamicList<label>& intersectedCells,
    const DynamicList<label>& submergedCells,
    const List<scalar>& shapeAlpha,
    Istream& is
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldDesc.type() != fieldType::typeName)
    {
        return false;
    }

    // Get value from stream
    const Type fieldValue = pTraits<Type>(is);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldDesc.name(),
        mesh.thisDb().time().timeName(),
        mesh.thisDb(),
        IOobject::MUST_READ
    );

    bool found(fieldHeader.typeHeaderOk<fieldType>(true));

    // Field exists
    if (found)
    {
        Info<< "    - set internal values of "
            << fieldHeader.headerClassName()
            << ": " << fieldDesc.name()
            << " = " << fieldValue << endl;

        fieldType field(fieldHeader, mesh, false);

        if (isNull(shapeAlpha))
        {
            field.primitiveFieldRef() = fieldValue;
        }
        else
        {
            if (intersectedCells.size())
            {
                forAll(intersectedCells, cellI)
                {
                    const Type oldValueI(field[intersectedCells[cellI]]);

                    field[intersectedCells[cellI]] =
                        shapeAlpha[cellI] * fieldValue
                      + (1.0 - shapeAlpha[cellI]) * oldValueI;
                }
            }

            if (submergedCells.size())
            {
                for (const label celli : submergedCells)
                {
                    field[celli] = fieldValue;
                }
            }
        }

        // Make boundary fields consistent - treat like zeroGradient
        for (auto& pfld : field.boundaryFieldRef())
        {
            pfld = pfld.patchInternalField();
        }

        // Handle any e.g. halo-swaps
        field.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << endl;
        }
    }
    else
    {
        Warning
            << "Field " << fieldDesc.name() << " not found" << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dispatcher for setting volume fields
struct setCellField
{
    autoPtr<setCellField> clone() const { return nullptr; }  // placeholder

    static bool apply
    (
        const fieldDescription& fieldDesc,
        const fvMesh& m,
        const DynamicList<label>& iCells,
        const DynamicList<label>& sCells,
        const List<scalar>& sAlpha,
        Istream& is
    )
    {
        return
        (
            setCellFieldType<scalar>(fieldDesc, m, iCells, sCells, sAlpha, is)
         || setCellFieldType<vector>(fieldDesc, m, iCells, sCells, sAlpha, is)
         || setCellFieldType<sphericalTensor>(fieldDesc, m, iCells, sCells, sAlpha, is)
         || setCellFieldType<symmTensor>(fieldDesc, m, iCells, sCells, sAlpha, is)
         || setCellFieldType<tensor>(fieldDesc, m, iCells, sCells, sAlpha, is)
        );
    }

    class iNew
    {
        const fvMesh& mesh_;
        const DynamicList<label>& intersectedCells_;
        const DynamicList<label>& submergedCells_;
        const List<scalar>& shapeAlpha_;

    public:

        iNew
        (
            const fvMesh& mesh,
            const DynamicList<label>& intersectedCells,
            const DynamicList<label>& submergedCells,
            const List<scalar>& shapeAlpha
        )
        :
            mesh_(mesh),
            intersectedCells_(intersectedCells),
            submergedCells_(submergedCells),
            shapeAlpha_(shapeAlpha)
        {}

        autoPtr<setCellField> operator()(Istream& is) const
        {
            const fieldDescription fieldDesc(is);

            bool ok = setCellField::apply
            (
                fieldDesc,
                mesh_,
                intersectedCells_,
                submergedCells_,
                shapeAlpha_,
                is
            );

            if (!ok)
            {
                ok = consumeUnused(fieldDesc, is);

                if (ok)
                {
                    // Not meant for us
                    Info<< "Skip " << fieldDesc.type()
                        << " for finite-volume" << nl;
                }
                else
                {
                    WarningInFunction
                        << "Unsupported field type: "
                        << fieldDesc.type() << endl;
                }
            }

            return nullptr;  // Irrelevant return value
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set field values using arbitrary shapes through a dictionary. For "
        "cells that are partially enclosed by these shapes, the values are "
        "averaged based on volume fractions."
    );

    // Disable and add some options
    argList::noParallel();          // Disable parallel function
    argList::noFunctionObjects();   // Don't use function objects

    argList::addOption("dict", "file", "Alternative setShapedFieldsDict");

    argList::addOption("np", "number", "Number of OpenMP threads");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "setOpenMP.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    scalar startTime(omp_get_wtime());

    #include "readControlDict.H"
    #include "setFieldValues.H"

    Info<< "Execution time: "
        << omp_get_wtime() - startTime
        << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
