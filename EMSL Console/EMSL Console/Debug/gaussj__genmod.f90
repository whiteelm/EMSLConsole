        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 25 17:49:38 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAUSSJ__genmod
          INTERFACE 
            SUBROUTINE GAUSSJ(N,ALPHA,BETA,B,T,W)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: B(N)
              REAL(KIND=8) :: T(N)
              REAL(KIND=8) :: W(N)
            END SUBROUTINE GAUSSJ
          END INTERFACE 
        END MODULE GAUSSJ__genmod
