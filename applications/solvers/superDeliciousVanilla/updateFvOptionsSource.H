    // Because we cannot directly call the fvc::flux operation on the
    // fvOptions(U) object, we have to extract the source term and store
    // it in a volVectorField.
    {
        fvVectorMatrix fvOptions_ = fvOptions(U);
        vectorField fvOptionsVelocitySource_ = fvOptions_.source();
        forAll(fvOptionsVelocitySource_, i)
        {
            fvOptionsVelocitySource[i] = fvOptionsVelocitySource_[i];
        }
        fvOptionsVelocitySource.correctBoundaryConditions();
    }
