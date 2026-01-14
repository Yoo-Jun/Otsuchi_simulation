#include "macro.h"

PROGRAM main
  USE misc
  USE geometry
  IMPLICIT NONE

  CALL init

  DO WHILE (t_current < t_end)
     IF (offline) THEN
        CALL step_offline
     ELSE
        CALL step
     END IF
  END DO

  CALL finalize

END PROGRAM main
