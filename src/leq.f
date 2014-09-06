       PROGRAM EQ
         IMPLICIT NONE
         
         INTEGER :: n
         REAL, DIMENSION(:, :), ALLOCATABLE :: A, B
         REAL, DIMENSION(:, :), ALLOCATABLE :: L, Lr, U, T, C
         REAL :: MIN_ERR
         REAL, DIMENSION(:, :), ALLOCATABLE :: X, X2, ERR
         
         INTERFACE
           FUNCTION MatrixProduct(A, B)
             REAL, INTENT(IN) :: A(:, :)
             REAL, INTENT(IN) :: B(:, :)

             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixProduct
           END FUNCTION MatrixProduct
           
           SUBROUTINE MatrixSeparate(A, L, U)
             REAL, INTENT(IN) :: A(:, :)
             REAL, INTENT(OUT) :: L(:, :)
             REAL, INTENT(OUT) :: U(:, :)
           END SUBROUTINE MatrixSeparate
           
           FUNCTION MatrixMinor(A, Row, Col)
             REAL, INTENT(IN) :: A(:, :)
             INTEGER, INTENT(IN) :: Row, Col
             
             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixMinor
           END FUNCTION MatrixMinor
           
           FUNCTION MatrixCofactor(A)
             REAL, INTENT(IN) :: A(:, :)
             
             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixCofactor
           END FUNCTION MatrixCofactor
           
           REAL FUNCTION MatrixDeterminant(A)
             REAL, INTENT(IN) :: A(:, :)
           END FUNCTION
           
           FUNCTION MatrixInverse(A)
             REAL, INTENT(IN) :: A(:, :)
             
             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixInverse
           END FUNCTION
         END INTERFACE
         
         WRITE (*,*), "Give dimension of Coefficients Matrix"
         READ (*,*), n
         
         ALLOCATE(A(n, n))
         ALLOCATE(B(n, 1))
         
         WRITE (*,*), "Give Coefficients Matrix (vertical): "
         READ (*, *), A
         WRITE (*,*), "Give right hand side Matrix: "
         READ (*, *), B
         
         CALL MatrixSeparate(A, L, U)
         
         ALLOCATE(Lr(n, n))
         ALLOCATE(T(n, n))
         ALLOCATE(C(n, 1))
         
         Lr = MatrixInverse(L)
         T = -MatrixProduct(Lr, U)
         C = MatrixProduct(Lr, B)
         
         ALLOCATE(X(n, 1))
         ALLOCATE(X2(n, 1))
         ALLOCATE(ERR(n, 1))
         
         X = 1
         MIN_ERR=1E-6
         ERR = MIN_ERR+1
         
         DO WHILE (MAXVAL(ERR)>MIN_ERR)
           X2 = MatrixProduct(T, X)+C
           ERR = Abs(X2-X)
           X = X2
           PRINT *, X
         END DO
         
         PRINT *, X
       END PROGRAM EQ
       
       RECURSIVE FUNCTION MatrixProduct(A, B) RESULT(Returning)
         REAL, INTENT(IN) :: A(:, :)
         REAL, INTENT(IN) :: B(:, :)
         
         REAL, DIMENSION(:, :), ALLOCATABLE :: Returning
         
         INTEGER :: DimA = Size(Shape(A)), DimB = Size(Shape(B))
         INTEGER, DIMENSION(:), ALLOCATABLE :: AShape, BShape
         REAL, DIMENSION(:, :), ALLOCATABLE :: FinalResult
         REAL, DIMENSION(:), ALLOCATABLE :: CurrRow
         REAL, DIMENSION(:), ALLOCATABLE :: CurrCol
         REAL, DIMENSION(:, :), ALLOCATABLE :: CurrRow2
         REAL, DIMENSION(:, :), ALLOCATABLE :: CurrCol2
         INTEGER :: I, J, Row, Col
         REAL :: TmpResult(1, 1)

         ALLOCATE(AShape(DimA))
         ALLOCATE(BShape(DimB))
         AShape = Shape(A)
         BShape = Shape(B)
         
         ALLOCATE(CurrRow(AShape(2)))
         ALLOCATE(CurrCol(BShape(1)))
         ALLOCATE(CurrRow2(1, AShape(2)))
         ALLOCATE(CurrCol2(BShape(1), 1))
         
         ALLOCATE(Returning(AShape(1), BShape(2)))
         ALLOCATE(FinalResult(AShape(1), BShape(2)))
         FinalResult=0
         
         IF (Size(Returning) == 1) THEN
           DO I=1,AShape(2)
             FinalResult(1,1) = FinalResult(1,1)+(A(1,I)*B(I,1))
           END DO
         ELSE
           DO Row=1,AShape(1)
             Do Col=1,BShape(2)
               CurrRow = A(Row, :)
               CurrCol = B(:, Col)
               I = AShape(2)
               J = BShape(1)
               CurrRow2 = Reshape(CurrRow, (/1, I/))
               CurrCol2 = Reshape(CurrCol,(/J, 1/))
               
               TmpResult = MatrixProduct(CurrRow2, CurrCol2)
               FinalResult(Row, Col) = TmpResult(1,1)
             END DO
           END DO
         END IF

         Returning = FinalResult
       END FUNCTION MatrixProduct
       
       SUBROUTINE MatrixSeparate(A, L, U)
         REAL, INTENT(IN) :: A(:, :)
         REAL, INTENT(OUT), ALLOCATABLE :: L(:, :)
         REAL, INTENT(OUT), ALLOCATABLE :: U(:, :)
           
         INTEGER, DIMENSION(2) :: AShape
         INTEGER :: Row, Col
         
         AShape=Shape(A)
           
         ALLOCATE(L(AShape(1), AShape(2)))
         ALLOCATE(U(AShape(1), AShape(2)))
           
         DO Row=1,AShape(1)
           DO Col=1,AShape(2)
             IF (Row<Col) THEN
               U(Row, Col) = A(Row, Col)
               L(Row, Col) = 0
             ELSE
               U(Row, Col) = 0
               L(Row, Col) = A(Row, Col)
             END IF
           END DO
         END DO
       END SUBROUTINE MatrixSeparate
       
       FUNCTION MatrixMinor(A, Row, Col)
         REAL, INTENT(IN) :: A(:, :)
         INTEGER, INTENT(IN) :: Row, Col
         
         REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixMinor
         REAL, DIMENSION(:, :), ALLOCATABLE :: FinalResult
         
         INTEGER, DIMENSION(2) :: AShape
         INTEGER :: iRow, iCol, fRow, fCol
         
         AShape = Shape(A)
         
         ALLOCATE(MatrixMinor(AShape(1)-1,AShape(2)-1))
         ALLOCATE(FinalResult(AShape(1)-1,AShape(2)-1))
         
         DO iRow=1,AShape(1)
           DO iCol=1,AShape(2)
             IF (iRow/=Row .AND. iCol/=Col) THEN
               IF (iRow<Row) THEN
                 fRow=iRow
               ELSE
                 fRow=iRow-1
               END IF
               
               IF (iCol<Col) THEN
                 fCol=iCol
               ELSE
                 fCol=iCol-1
               END IF
               
               FinalResult(fRow, fCol) = A(iRow, iCol)
             END IF
           END DO
         END DO
         
         MatrixMinor = FinalResult
       END FUNCTION MatrixMinor
       
       FUNCTION MatrixCofactor(A)
         REAL, INTENT(IN) :: A(:, :)

         REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixCofactor
         
         REAL, DIMENSION(:, :), ALLOCATABLE :: FinalResult
         REAL, DIMENSION(:, :), ALLOCATABLE :: M
         REAL :: D
         
         INTEGER, DIMENSION(2) :: AShape
         INTEGER :: iRow, iCol
         
         INTERFACE
           FUNCTION MatrixMinor(A, Row, Col)
             REAL, INTENT(IN) :: A(:, :)
             INTEGER, INTENT(IN) :: Row, Col
             
             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixMinor
           END FUNCTION MatrixMinor
           
           REAL FUNCTION MatrixDeterminant(A)
             REAL, INTENT(IN) :: A(:, :)
           END FUNCTION
         END INTERFACE
         
         AShape = Shape(A)
         
         ALLOCATE(MatrixCofactor(AShape(1),AShape(2)))
         ALLOCATE(FinalResult(AShape(1),AShape(2)))
         
         IF (AShape(1)==1 .AND. AShape(2)==1) THEN
           FinalResult(1, 1) = 1
         ELSE
           ALLOCATE(M(AShape(1)-1,AShape(2)-1))
           
           DO iRow=1,AShape(1)
             DO iCol=1,AShape(2)
               M = MatrixMinor(A, iRow, iCol)
               D = MatrixDeterminant(M)
               FinalResult(iRow, iCol) = ((-1)**(iRow+iCol))*D
             END DO
           END DO
         END IF
         
         MatrixCofactor = FinalResult
       END FUNCTION MatrixCofactor
       
       REAL FUNCTION MatrixDeterminant(A)
         REAL, INTENT(IN) :: A(:, :)
         
         REAL :: FinalResult
         REAL, DIMENSION(:, :), ALLOCATABLE :: C
         
         INTEGER, DIMENSION(2) :: AShape
         INTEGER :: iRow
         
         INTERFACE
           FUNCTION MatrixCofactor(A)
             REAL, INTENT(IN) :: A(:, :)
             
             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixCofactor
           END FUNCTION MatrixCofactor
         END INTERFACE
         
         AShape = Shape(A)
         FinalResult=0
         
         IF (AShape(1)==1 .AND. AShape(2)==1) THEN
           FinalResult = A(1, 1)
         ELSE
           ALLOCATE(C(AShape(1),AShape(2)))
           C=MatrixCofactor(A)
           
           DO iRow=1,AShape(1)
             FinalResult = FinalResult+(A(iRow, 1)*C(iRow, 1))
           END DO
         END IF
         
         MatrixDeterminant = FinalResult
       END FUNCTION MatrixDeterminant
       
       FUNCTION MatrixInverse(A)
         REAL, INTENT(IN) :: A(:, :)
         
         REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixInverse
         REAL, DIMENSION(:, :), ALLOCATABLE :: FinalResult
         
         REAL, DIMENSION(:, :), ALLOCATABLE :: C
         REAL, DIMENSION(:, :), ALLOCATABLE :: CT
         REAL :: D
         
         INTEGER, DIMENSION(2) :: AShape
         
         INTERFACE
           FUNCTION MatrixCofactor(A)
             REAL, INTENT(IN) :: A(:, :)
             
             REAL, DIMENSION(:, :), ALLOCATABLE :: MatrixCofactor
           END FUNCTION MatrixCofactor
           
           REAL FUNCTION MatrixDeterminant(A)
             REAL, INTENT(IN) :: A(:, :)
           END FUNCTION
         END INTERFACE
         
         AShape = Shape(A)
         
         ALLOCATE(MatrixInverse(AShape(1),AShape(2)))
         ALLOCATE(FinalResult(AShape(1),AShape(2)))
         
         ALLOCATE(C(AShape(1),AShape(2)))
         ALLOCATE(CT(AShape(1),AShape(2)))
         C = MatrixCofactor(A)
         CT = Transpose(C)
         D = MatrixDeterminant(A)
         
         FinalResult = (1/D)*CT
         
         MatrixInverse = FinalResult
       END FUNCTION