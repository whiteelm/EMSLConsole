        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 25 17:49:38 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SCOUTP__genmod
          INTERFACE 
            SUBROUTINE SCOUTP(N,C,Z,WC,W,BETAM,NPTSQ)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: C
              COMPLEX(KIND=8) :: Z(N)
              COMPLEX(KIND=8) :: WC
              COMPLEX(KIND=8) :: W(N)
              REAL(KIND=8) :: BETAM(N)
              INTEGER(KIND=4) :: NPTSQ
            END SUBROUTINE SCOUTP
          END INTERFACE 
        END MODULE SCOUTP__genmod
