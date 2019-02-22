﻿using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: perhaps the solver should expose the assembler, instead of wrapping it. The assembler's interface would have to be 
//      simplified a bit though. That would violate LoD, but so does MSolve in general.
namespace ISAAR.MSolve.Solvers
{
    /// <summary>
    /// Helps translate the physical model into a linear system and then solves the latter. 
    /// </summary>
    public interface ISolver_v2 : ISystemMatrixObserver
    {
        /// <summary>
        /// A dictionary that maps subdomain ids to linear systems.
        /// </summary>
        IReadOnlyDictionary<int, ILinearSystem_v2> LinearSystems { get; }

        /// <summary>
        /// Assembles the matrix that corresponds to the free freedom degrees of the whole subdomain from the matrices of its 
        /// elements.
        /// </summary>
        /// <param name="subdomain">The subdomain whose corresponding matrix will be assembled.</param>
        /// <param name="elementMatrixProvider">
        /// Determines the matrix calculated for each element (e.g. stiffness, mass, etc.)
        /// </param>
        IMatrix BuildGlobalMatrix(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider); //TODO: Ideally the provider/analyzer will not even have to pass the subdomain.

        /// <summary>
        /// Assembles the matrices that correspond to the free and constrained freedom degrees of the whole subdomain 
        /// from the matrices of its elements. If we denote the matrix as A, the free dofs as f and the constrained dofs as c
        /// then: A = [ Aff Acf^T; Acf Acc ] (Matlab notation). This method returns Aff, Acf, Acc.
        /// </summary>
        /// <param name="subdomain">The subdomain whose corresponding matrix will be assembled.</param>
        /// <param name="elementMatrixProvider">
        /// Determines the matrix calculated for each element (e.g. stiffness, mass, etc.)
        /// </param>
        (IMatrix matrixFreeFree, IMatrix matrixConstrFree, IMatrix matrixConstrConstr) 
            BuildGlobalSubmatrices(ISubdomain_v2 subdomain, IElementMatrixProvider_v2 elementMatrixProvider);

        /// <summary>
        /// Initializes the state of this <see cref="ISolver_v2"/> instance. This needs to be called only once, since it  
        /// could potentially perform actions that must not be repeated or are too expensive.
        /// </summary>
        void Initialize();

        /// <summary>
        /// Solves multiple linear systems A * X = B, where: A is one of the matrices stored in <see cref="LinearSystems"/>,
        /// B is the corresponding matrix in <paramref name="otherMatrix"/> and X is the corresponding matrix that will be 
        /// calculated as the result of inv(A) * B. 
        /// </summary>
        /// <param name="otherMatrix">
        /// The right hand side matrix for each subdomain. If the linear systems are A * X = B, then B is one of the matrices in
        /// <paramref name="otherMatrix"/>.</param>
        Dictionary<int, Matrix> InverseSystemMatrixTimesOtherMatrix(Dictionary<int, IMatrixView> otherMatrix);

        /// <summary>
        /// Orders the free and optionally the constrained freedom degrees of the model. Also remember to reset the linear 
        /// systems. 
        /// </summary>
        /// <param name="alsoOrderConstrainedDofs">
        /// If true, the constrained dofs will also be ordered. Else, <see cref="ISubdomain_v2.ConstrainedDofOrdering"/> will
        /// be null.
        /// </param>
        void OrderDofs(bool alsoOrderConstrainedDofs);
        //TODO: Would it be better if the solver didn't modify the model? It would return the ordering 
        //      and the analyzer would. However the assembler should still be notified.
        //TODO: I think this should also reset the linear systems. Is there any chance that OrderDofs() should be called but 
        //      ILinearSystem.Reset() shouldn't?

        /// <summary>
        /// Notifies this <see cref="ISolver_v2"/> that it cannot overwrite the data of <see cref="ILinearSystem_v2.Matrix"/>.
        /// Some solvers would otherwise overwrite the matrices (e.g. with the factorization) to avoid using extra memory.
        /// </summary>
        void PreventFromOverwrittingSystemMatrices();

        /// <summary>
        /// Solves the linear systems.
        /// </summary>
        void Solve();
    }
}
