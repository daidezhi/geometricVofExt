// // Set default values
// ITstream& is(defaultValues.stream());

// is.fatalCheck(FUNCTION_NAME);

// token tok(is);

// is >> tok;

// while (!tok.isPunctuation(token::END_LIST))
// {
//     is.putBack(tok);

//     word fieldTypeName(tok.wordToken());

//     if (fieldTypeName.ends_with("Value"))
//     {
//         fieldTypeName.erase(fieldTypeName.size()-5);
//     }

//     if (fieldTypeName == volScalarField::typeName)
//     {
//         const fieldDescription<scalar> fieldDesc(is);
//         Info << fieldDesc.name() << " ";
//         Info << fieldDesc.type() << " ";
//         Info << fieldDesc.value() << endl;
//     }
//     else if (fieldTypeName == volVectorField::typeName)
//     {
//         const fieldDescription<vector> fieldDesc(is);
//         Info << fieldDesc.name() << " ";
//         Info << fieldDesc.type() << " ";
//         Info << fieldDesc.value() << endl;
//     }

//     is >> tok;
// }

// Set default values
{
    const entry* eptr
    (
        setShapedFieldsDict.findEntry("defaultFieldValues", keyType::LITERAL)
    );

    if (eptr)
    {
        ITstream& is(eptr->stream());

        Info << "Setting volume field default values" << endl;

        PtrList<setCellField> defaultFieldValues
        (
            is,
            setCellField::iNew
            (
                mesh,
                DynamicList<label>(),
                DynamicList<label>(),
                List<scalar>::null()
            )
        );

        Info << endl;
    }
}


// Set shaped values
{
    const scalarField& cellVolumes(mesh.cellVolumes());

    for (const entry& shapei: shapes)
    {
        const word shapeIName(shapei.dict().dictName());
        const word shapeIType(shapei.dict().get<word>("type"));

        Info<< "Setting field values with shape: " << shapeIName
            << " (type: " << shapeIType << ")"
            << endl;

        const geometricVofExt::CGALVof::CGALSurfaceMesh CGALSMeshI
        (
            geometricVofExt::CGALVof::CGALSurfaceMesh(shapei.dict())
        );

        if (!CGALSMeshI.surface().is_empty())
        {
            Info << "    - calculate mesh point and cell locations" << endl;

            const geometricVofExt::CGALVof::cellToSurfaceMeshLocation cellLocs
            (
                mesh,
                CGALSMeshI,
                extraBooleanCheck
            );

            const DynamicList<label>& intersectedCells(cellLocs.intersectedCells());
            const DynamicList<label>& submergedCells(cellLocs.submergedCells());
            //const DynamicList<label>& aridCells(cellLocs.aridCells());

            Info << "    - calculating fraction values" << endl;

            const label nIntersectedCells(intersectedCells.size());
            const label nSubmergedCells(submergedCells.size());

            List<scalar> shapeAlpha(nIntersectedCells);
            
            #pragma omp parallel for
            forAll(intersectedCells, cellI)
            {
                const geometricVofExt::CGALVof::CGALSurfaceMesh cCell
                (
                    mesh,
                    intersectedCells[cellI],
                    triangulationScheme
                );

                geometricVofExt::CGALVof::surfaceMeshBooleanOperation booleanOp
                (
                    CGALSMeshI,
                    cCell,
                    word("intersection")
                );

                booleanOp.compBooleanOp();

                shapeAlpha[cellI] =
                    booleanOp.volume() / cellVolumes[intersectedCells[cellI]];
                
                shapeAlpha[cellI] = min(1.0, max(0.0, shapeAlpha[cellI]));
            }

            const entry* eptr
            (
                shapei.dict().findEntry("fieldValues", keyType::LITERAL)
            );

            if (eptr)
            {
                ITstream& is(eptr->stream());

                PtrList<setCellField> shapedFieldValues
                (
                    is,
                    setCellField::iNew
                    (
                        mesh,
                        intersectedCells,
                        submergedCells,
                        shapeAlpha
                    )
                );
            }
        }

        Info << endl;
    }
}
