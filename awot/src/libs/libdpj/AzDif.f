C  "<08-Oct-1993 21:54:25><UTC>"
*************************************************************************
      Real*4 Function AzDif(arg1,arg2, difMax)

*  Input:
*     arg1 and arg2 should be between 0.0 and 2.*difMax.
*  Output:
*     case 1: if  -difMax <= arg1 -arg2 <= difMax,
*                then (arg1-arg2) will be returned.
*     case 2: if   arg1 -arg2 < -difMax,
*                then (arg1-arg2)+2*difMax will be returned.
*     case 3: if   arg1 -arg2 > difMax,
*                then (arg1-arg2)-2*difMax will be returned.
*  I.e., return difference of two arguments, where modulo
*  arithmetic of 2*difMax means the difference has magnitude <= difMax.

      Real*4 arg1, arg2, difMax
      Real*4 dif
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      dif = arg1 - arg2
      if (dif  .lt. -difMax)then
	  AzDif = dif + 2.*difMax
	  return
      else if (dif .gt. difMax)then
	  AzDif = dif - 2.*difMax
	  return
      end if
      AzDif = dif
      return 
      end ! Real*4 function azDif ends 
