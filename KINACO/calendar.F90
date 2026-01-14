#include "macro.h"

MODULE calendar
  USE misc
  IMPLICIT NONE
  PRIVATE
  PUBLIC dtime, idtime, idtime_r4
  PUBLIC datetime_seconds, interval_seconds, format_datetime, day_per_year, strip_date, strip_time
  PUBLIC t_current, t_start, t_end, current_datetime, start_datetime, end_datetime
  PUBLIC seconds2date, seconds2time, doy2dom
  PUBLIC calendar_type, calendar_format

  INTEGER, SAVE :: calendar_type   = 0  ! 0: 360 days per year
                                        ! 1: 365 days per year
                                        ! 2: Gregorian

  INTEGER, SAVE :: calendar_format = 0  ! 0: YYYYDDD_HHMMSS   (DDD: day_of_year)
                                        ! 1: YYYYmmdd_HHMMSS  (mm: month,      dd: day_of_month)
                                        ! 2: YYYYmmmdd_HHMMSS (mmm:month_abbr, dd:day_of_month)

  CHARACTER(3), PARAMETER ::    abbr_month(12) = (/'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
  INTEGER,      PARAMETER :: day_per_month(12) = (/   31,    28,    31,    30,    31,    30,    31,    31,    30,    31,    30,    31/)

  REAL(4), SAVE :: dtime  = 0.0
  REAL(8), SAVE :: idtime = 0.0
  REAL(4), SAVE :: idtime_r4 = 0.0

  REAL(8), SAVE :: t_current
  REAL(8), SAVE :: t_start
  REAL(8), SAVE :: t_end

  CHARACTER(16), SAVE :: current_datetime, start_datetime, end_datetime

CONTAINS
  REAL(8) FUNCTION datetime_seconds(datetime)
    CHARACTER(*), INTENT(IN) :: datetime
    INTEGER :: fmt
    INTEGER :: year, doy, hour, minute, second
    INTEGER :: month, dom
    INTEGER :: dpm(0:12)
    INTEGER :: i, m, n
    INTEGER :: iostat

    CHARACTER(16), SAVE :: previous_arg    = ""
    REAL(8),       SAVE :: previous_result = UNDEF

    IF (trim(datetime) == trim(previous_arg) .AND. previous_result /= UNDEF) THEN
       datetime_seconds = previous_result
       RETURN
    END IF

    previous_arg = trim(datetime)

    fmt = inquire_format(datetime)

    SELECT CASE (fmt)
    CASE (-2)
       CALL assert(.FALSE., "'*S/M/H/D'-style DATETIME format is valid only for INTERVAL")
    CASE (-1)
       m = 0
    CASE(0)
       READ(datetime(1:4), '(I4.4)') year
       READ(datetime(5:7), '(I3.3)') doy
       m = 8
    CASE (1)
       READ(datetime(1:4), '(I4.4)') year
       READ(datetime(5:6), '(I2.2)') month
       READ(datetime(7:8), '(I2.2)') dom
       m = 9
    CASE (2)
       READ(datetime(1:4), '(I4.4)') year
       DO month=1, 12
          IF (datetime(5:7) == abbr_month(month)) EXIT
       END DO
       READ(datetime(8:9), '(I2.2)') dom
       m = 10
    END SELECT

    IF (fmt==1 .OR. fmt==2) THEN
       dpm(0) = 0
       IF (calendar_type==0) THEN
          dpm(1:12) = 30
       ELSE
          dpm(1:12) = day_per_month(:)
          IF (calendar_type==2 .AND. is_leap_year(year)) dpm(2) = 29
       END IF
       CALL assert(month >=1 .AND. month <= 12,      "invalid month in DATETIME '"//trim(datetime)//"'")
       CALL assert(dom >= 1 .AND. dom <= dpm(month), "invalid day of month in DATETIME '"//trim(datetime)//"'")
       doy = dom + sum(dpm(0:month-1)) - 1
    END IF

    IF (len_trim(datetime) == m-1) THEN
       hour   = 0
       minute = 0
       second = 0
    ELSE
       READ(datetime(m+1:m+2), '(I2.2)') hour
       READ(datetime(m+3:m+4), '(I2.2)') minute
       READ(datetime(m+5:m+6), '(I2.2)') second
    END IF

    CALL assert(second >= 0 .AND. second < 60, "invalid second in DATETIME '"// trim(datetime)//"'")
    CALL assert(minute >= 0 .AND. minute < 60, "invalid minute in DATETIME '"// trim(datetime)//"'")
    CALL assert(hour   >= 0 .AND. hour   < 24, "invalid hour in DATETIME '"// trim(datetime)//"'")

    CALL assert(doy >= 0 .AND. doy < day_per_year(year), "invalid doy in DATETIME '"//trim(datetime)//"'")

    CALL assert(year >= 0 .AND. year < 10000, "invalid year in DATETIME '"//trim(datetime)//"'")
    IF (calendar_type==2) THEN
       datetime_seconds = DBLE(year*365 + count_leap_year(year-1)) * 24*60*60
    ELSE
       datetime_seconds = DBLE(year) * day_per_year(0)*24*60*60
    END IF

    datetime_seconds = datetime_seconds + DBLE(doy)    * 24*60*60
    datetime_seconds = datetime_seconds + DBLE(hour)   * 60*60
    datetime_seconds = datetime_seconds + DBLE(minute) * 60
    datetime_seconds = datetime_seconds + DBLE(second)

    previous_result = datetime_seconds

  END FUNCTION datetime_seconds

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) FUNCTION interval_seconds(interval)
    CHARACTER(*), INTENT(IN) :: interval

    INTEGER    :: l, i
    CHARACTER  :: c
    INTEGER(8) :: n

    SELECT CASE (inquire_format(interval))
    CASE (-2)
       l = len_trim(interval)
       READ(interval(1:l-1), *) n
       SELECT CASE (interval(l:l))
       CASE ('S', 's')
          interval_seconds = REAL(n,8)
       CASE ('M', 'm')
          interval_seconds = REAL(n,8) * 60
       CASE ('H', 'h')
          interval_seconds = REAL(n,8) * 3600
       CASE ('D', 'd')
          interval_seconds = REAL(n,8) * 86400
       CASE ('Y', 'y')
          interval_seconds = REAL(n,8) * 86400 * day_per_year(0)
       CASE ('T', 't')
          interval_seconds = REAL(n,8) * dtime
       CASE DEFAULT
          CALL assert(.FALSE., "invalid INTERVAL '"//trim(interval)//"'")
       END SELECT
    CASE (-1,0)
       interval_seconds = datetime_seconds(interval)
    CASE DEFAULT
          CALL assert(.FALSE., "invalid INTERVAL '"//trim(interval)//"'")
    END SELECT

    CALL assert(interval_seconds >= 0, "invalid INTERVAL '"//trim(interval)//"'")

  END FUNCTION interval_seconds

!-----------------------------------------------------------------------------------------------------------------------

  PURE SUBROUTINE seconds2date(seconds, year, doy, rem)
    REAL(8), INTENT(IN)  :: seconds
    INTEGER, INTENT(OUT) :: year
    INTEGER, INTENT(OUT) :: doy
    REAL(8), INTENT(OUT), OPTIONAL :: rem

    INTEGER :: n

    doy  = INT(seconds / 86400)
    IF (present(rem)) rem = seconds - doy*86400D0

    SELECT CASE (calendar_type)
    CASE (0)
       year = INT(doy / 360)
       doy  = doy - year*360

    CASE (1)
       year = INT(doy / 365)
       doy  = doy - year*365

    CASE(2)
       year=0

       IF (doy >= 365) THEN
          doy  = doy - 365 !adjust for unused year-zero
          year = year + 1

          n = int(doy / (365*400+97))
          doy  = doy  - n*(365*400+97)
          year = year + n*400

          n = min(int(doy / (365*100+24)), 3)
          doy  = doy  - n*(365*100+24)
          year = year + n*100

          n    = int(doy / (365*4+1))
          doy  = doy  - n*(365*4+1)
          year = year + n*4

          n    = min(int(doy / 365), 3)
          doy  = doy  - n*365
          year = year + n
       END IF
    CASE DEFAULT
       year = -9999
       doy  = 0
    END SELECT

  END SUBROUTINE seconds2date

  PURE SUBROUTINE seconds2time(seconds, hour, min, sec, rem)
    REAL(8), INTENT(IN)  :: seconds
    INTEGER, INTENT(OUT) :: hour
    INTEGER, INTENT(OUT) :: min
    INTEGER, INTENT(OUT) :: sec
    REAL(8), INTENT(OUT), OPTIONAL :: rem

    REAL(8) :: tmp

    tmp = seconds - INT(seconds/86400)*86400D0

    hour = INT(tmp / 3600)
    tmp  = tmp - hour*3600D0
    min  = INT(tmp / 60)
    tmp  = tmp - min*60D0

    sec = INT(tmp)

    IF (present(rem)) rem = tmp - sec

  END SUBROUTINE seconds2time


!-----------------------------------------------------------------------------------------------------------------------


  CHARACTER(20) FUNCTION format_datetime(t, fmt, omit_date, omit_time, fraction)
    REAL(8), INTENT(IN) :: t
    INTEGER, INTENT(IN), OPTIONAL :: fmt
    LOGICAL, INTENT(IN), OPTIONAL :: omit_time, omit_date, fraction
    INTEGER, SAVE :: year = -1
    INTEGER, SAVE :: doy  = -1
    INTEGER, SAVE :: hour = -1
    INTEGER, SAVE :: min  = -1
    INTEGER, SAVE :: sec  = -1
    REAL(8), SAVE :: previous_t = UNDEF

    INTEGER :: fmt_
    LOGICAL :: ot_, od_
    CHARACTER(3) :: frac_

    IF (t /= previous_t) THEN
       previous_t = t

       CALL seconds2date(t, year, doy)
       CALL seconds2time(t, hour, min, sec)
    END IF

    fmt_ = calendar_format
    IF (present(fmt)) fmt_ = fmt

    od_ = .FALSE.
    IF (present(omit_date)) od_ = omit_date

    IF (od_) THEN
       format_datetime = ""
    ELSE
       format_datetime = trim(format_date(year, doy, fmt_))
    END IF

    ot_ = .FALSE.
    IF (present(omit_time)) ot_ = omit_time

    IF (.NOT. ot_) format_datetime = trim(format_datetime) // '_' // format_time(hour, min, sec)

    IF (present(fraction)) THEN
       IF (fraction) THEN
          WRITE(frac_,'(I3.3)') INT((t - FLOOR(t))/dtime)
          format_datetime = trim(format_datetime) // '_' // frac_
       END IF
    END IF


  END FUNCTION format_datetime

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(9) FUNCTION strip_date(datetime)
    CHARACTER(*), INTENT(IN) :: datetime

    INTEGER :: fmt

    fmt = inquire_format(datetime)

    SELECT CASE (fmt)
    CASE (0)
       IF (len_trim(datetime)==7) THEN
          strip_date = trim(datetime)
       ELSE
          strip_date = datetime(1:7)
       END IF
    CASE (1)
       IF (len_trim(datetime)==8) THEN
          strip_date = trim(datetime)
       ELSE
          strip_date = datetime(1:8)
       END IF
    CASE (2)
       IF (len_trim(datetime)==9) THEN
          strip_date = trim(datetime)
       ELSE
          strip_date = datetime(1:9)
       END IF
    CASE DEFAULT
       CALL assert(.FALSE., "'"//trim(datetime)//"'"//" is not valid datetime format")
    END SELECT

  END FUNCTION strip_date

  CHARACTER(6) FUNCTION strip_time(datetime)
    CHARACTER(*), INTENT(IN) :: datetime

    INTEGER :: fmt

    fmt = inquire_format(datetime)

    SELECT CASE (fmt)
    CASE (-1)
       strip_time = trim(datetime)
    CASE (0)
       IF (len_trim(datetime)==7) THEN
          strip_time = "000000"
       ELSE
          strip_time = datetime(9:14)
       END IF
    CASE (1)
       IF (len_trim(datetime)==8) THEN
          strip_time = "000000"
       ELSE
          strip_time = datetime(10:15)
       END IF
    CASE (2)
       IF (len_trim(datetime)==9) THEN
          strip_time = "000000"
       ELSE
          strip_time = datetime(11:16)
       END IF
    CASE DEFAULT
       CALL assert(.FALSE., "'"//trim(datetime)//"'"//" is not valid datetime format")
    END SELECT

  END FUNCTION strip_time

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION day_per_year(year)
    INTEGER, INTENT(IN) :: year

    SELECT CASE (calendar_type)
    CASE (0)
       day_per_year = 360
    CASE (1)
       day_per_year = 365
    CASE (2)
       IF (is_leap_year(year)) THEN
          day_per_year = 366
       ELSE
          day_per_year = 365
       END IF
    END SELECT
  END FUNCTION day_per_year

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION is_leap_year(year)
    INTEGER, INTENT(IN) :: year

    is_leap_year = (year > 0) .AND. (mod(year,4)==0) .AND. (mod(year,100)/=0 .OR. mod(year,400)==0)

  END FUNCTION is_leap_year

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION count_leap_year(year)
    INTEGER, INTENT(IN) :: year

    count_leap_year = max(0, int(year/4) - int(year/100) + int(year/400))

  END FUNCTION count_leap_year

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER FUNCTION inquire_format(datetime)
    CHARACTER(*), INTENT(IN) :: datetime

    INTEGER :: l, i
    LOGICAL :: flag

    l = len_trim(datetime)

    CALL assert(l > 0, "empty string passed to CALENDAR%INQUIRE_FORMAT")

    IF (datetime(l:l) == 's' .OR. datetime(l:l) == 'S' .OR. &
        datetime(l:l) == 'm' .OR. datetime(l:l) == 'M' .OR. &
        datetime(l:l) == 'h' .OR. datetime(l:l) == 'H' .OR. &
        datetime(l:l) == 'd' .OR. datetime(l:l) == 'D' .OR. &
        datetime(l:l) == 'y' .OR. datetime(l:l) == 'Y' .OR. &
        datetime(l:l) == 't' .OR. datetime(l:l) == 'T') THEN

       CALL assert(is_number(trim(adjustl(datetime(1:l-1)))), "invalid DATETIME '"//trim(datetime)//"'")
       inquire_format = -2
    ELSE
       SELECT CASE (l)
       CASE(6)
          CALL assert(is_number(datetime(1:6)), "invalid DATETIME '"//trim(datetime)//"'")
          inquire_format = -1
       CASE(7,14)
          CALL assert(is_number(datetime(1:7)), "invalid DATETIME '"//trim(datetime)//"'")
          IF (l==14) CALL assert(datetime(8:8)=='_' .AND. is_number(datetime(9:14)),  "invalid DATETIME '"//trim(datetime)//"'")
          inquire_format = 0
       CASE(8,15)
          CALL assert(is_number(datetime(1:8)), "invalid DATETIME '"//trim(datetime)//"'")
          IF (l==15) CALL assert(datetime(9:9)=='_' .AND. is_number(datetime(10:15)), "invalid DATETIME '"//trim(datetime)//"'")
          inquire_format = 1
       CASE(9,16)
          CALL assert (is_number(datetime(1:4)) .AND.is_number(datetime(8:9)), "invalid DATETIME '"//trim(datetime)//"'")
          flag = .FALSE.
          DO i=1, 12
             flag = datetime(5:7)==abbr_month(i)
             IF (flag) EXIT
          END DO
          CALL assert(flag, "invalid DATETIME '"//trim(datetime)//"'")
          IF (l==16) CALL assert(datetime(10:10)=='_' .AND. is_number(datetime(11:16)), "invalid DATETIME '"//trim(datetime)//"'")
          inquire_format = 2
       CASE DEFAULT
          CALL assert(.FALSE., "invalid DATETIME '"//trim(datetime)//"'")
          inquire_format = -99
       END SELECT
    END IF

  END FUNCTION inquire_format

!-----------------------------------------------------------------------------------------------------------------------

  PURE SUBROUTINE doy2dom(year, doy, mon, dom)
    INTEGER, INTENT(IN)  :: year, doy
    INTEGER, INTENT(OUT) :: mon,  dom

    INTEGER :: dpm(0:12)

    dpm(0) = 0
    IF (calendar_type == 0) THEN
       dpm(1:12) = 30
    ELSE
       dpm(1:12) = day_per_month(:)
       IF (calendar_type == 2 .AND. is_leap_year(year)) dpm(2) = 29
    END IF

    DO mon = 1, 12
       IF (doy+1 <= sum(dpm(0:mon))) THEN
          dom = doy+1 - sum(dpm(0:mon-1))
          EXIT
       END IF
    END DO

    !CALL assert(sum(dpm(:)) == day_per_year(year), "internal error in CALENDAR%DOY2DOM")
    !CALL assert(mon <= 12,                         "internal error in CALENDAR%DOY2DOM")
    !CALL assert(dom >= 1 .AND. dom <= dpm(mon),    "internal error in CALENDAR%DOY2DOM")

  END SUBROUTINE doy2dom

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(9) FUNCTION format_date(year, doy, fmt)
    INTEGER, INTENT(IN) :: year, doy
    INTEGER, INTENT(IN) :: fmt

    INTEGER :: mon, dom

    CHARACTER(9) :: result

    CALL assert(year >= 0 .AND. year < 100000,             "invalid YEAR")
    CALL assert(doy  >= 0 .AND. doy  < day_per_year(year), "invalid DOY")

    result = ""

    SELECT CASE (fmt)
    CASE (0)
    CASE (1:2)
       CALL doy2dom(year, doy, mon, dom)
    CASE DEFAULT
       CALL assert(.FALSE., "unssupported calndar format")
    END SELECT

    WRITE(result(1:4), '(I4.4)') year
    SELECT CASE (fmt)
    CASE (0)
       WRITE(result(5:7), '(I3.3)') doy
    CASE (1)
       WRITE(result(5:6), '(I2.2)') mon
       WRITE(result(7:8), '(I2.2)') dom
    CASE (2)
       result(5:7) = abbr_month(mon)
       WRITE(result(8:9), '(I2.2)') dom
    END SELECT

    format_date = trim(result)

  END FUNCTION format_date

  CHARACTER(6) FUNCTION format_time(hour, min, sec)
    INTEGER, INTENT(IN) :: hour, min, sec

    CHARACTER(6) :: result

    CALL assert(hour >= 0 .AND. hour < 24, "invalid HOUR")
    CALL assert(min  >= 0 .AND. min  < 60, "invalid MIN")
    CALL assert(sec  >= 0 .AND. sec  < 60, "invalid SEC")

    WRITE(result(1:2), '(I2.2)') hour
    WRITE(result(3:4), '(I2.2)') min
    WRITE(result(5:6), '(I2.2)') sec

    format_time = result

  END FUNCTION format_time

END MODULE calendar
