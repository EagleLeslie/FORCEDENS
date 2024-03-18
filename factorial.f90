RECURSIVE FUNCTION factorial(n) RESULT (fac_result)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: fac_result

  IF (n <= 1) THEN
    fac_result = 1
  ELSE
    fac_result = n * factorial(n-1)
  END IF

END FUNCTION factorial
