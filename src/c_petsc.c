#include <complex.h>

#include <petscksp.h>

int petsc_solve(int n_, double complex *A_, double complex *b_, double complex *x_)
{
    Vec            x, b;
    Mat            A;
    KSP            ksp;
    PC             pc;
    PetscReal      norm, tol=1.e-14;
    PetscErrorCode ierr;
    PetscInt       i, j, n = n_, col[n_], its;
    PetscScalar    neg_one = -1.0, one = 1.0, value[n_], *x_array;

    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    printf("  Converting matrix to PETSc\n");
    for (i=0; i<n; i++) col[i] = i;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) value[j] = A_[i*n+j];
        ierr = MatSetValues(A,1,&i,n,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    printf("    Done\n");
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    for (i=0; i<n; i++) value[i] = b_[i];
    ierr = VecSetValues(b, n, col, value, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

    ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    // For a full list, see:
    // http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html
    ierr = KSPSetType(ksp, KSPCGS);CHKERRQ(ierr);

    printf("  Solving...\n");
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    printf("  Done\n");

//    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    ierr = VecGetArray(x, &x_array);CHKERRQ(ierr);
    for (i=0; i<n; i++) x_[i] = x_array[i];

    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Iterations %D\n", its);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

    return 0;
}

int petsc_init()
{
    char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

    PetscErrorCode ierr;
    PetscMPIInt    size;
    int argc = 0;
    char **args;

    PetscInitialize(&argc,&args,(char *)0,help);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(PETSC_COMM_WORLD, 1, "This only works on one processor.");
}

int petsc_finalize()
{
    PetscErrorCode ierr;
    ierr = PetscFinalize();
}
