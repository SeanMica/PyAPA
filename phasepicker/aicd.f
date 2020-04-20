C FILE: AICD.F
       SUBROUTINE aicd(aic,aicder,data1,npts)
C
C      Calculate AIC for picker
C
       INTEGER npts
       REAL*8 data1(npts)
       REAL*8 aic(npts)
       REAL*8 aicder(npts)
       REAL*8 data2(npts)
       
       REAL*8 sum_data, sum_data2
       REAL*8 sum_data_k, mean_data_k1, mean_data_k2
       REAL*8 sum_data2_k, mean_data2_k1, mean_data2_k2
       
       
Cf2py  intent(in) npts
Cf2py  intent(in) data1
Cf2py  intent(out) aic
Cf2py  intent(out) aicder
Cf2py  depend(npts) aic
Cf2py  depend(npts) aicder
Cf2py  depend(npts) data1

       data2 = 0.d0
       sum_data = 0.d0
       sum_data2 = 0.d0

       DO i = 1, npts
          data2(i)  = data1(i)*data1(i)
          sum_data  = sum_data + data1(i)
          sum_data2 = sum_data2 + data2(i)
       ENDDO
       
       sum_data_k = sum_data - data1(npts)
       sum_data2_k = sum_data2 - data2(npts)
       x=0.0
       
       DO i = npts - 1, 2, -1
          sum_data_k = sum_data_k - data1(i)
          mean_data_k1 = (sum_data_k / i)*(sum_data_k / i)
          mean_data_k2 = ((sum_data-sum_data_k)/(npts-i))**2
          
          sum_data2_k = sum_data2_k - data2(i)
          mean_data2_k1 = sum_data2_k/i
          mean_data2_k2 = (sum_data2 - sum_data2_k)/(npts-i)
          
          a = i*log(mean_data2_k1-mean_data_k1)+(npts-i)*
     &          log(mean_data2_k2-mean_data_k2)
          
          IF ( a == -log(x) ) THEN
            a = aic(i+1)
          ENDIF
          
          aic(i) = a
          
       ENDDO
       
       aic(1) = aic(2)
       aic(npts) = aic(npts-1)
       
       DO i = 2, npts
          aicder(i) = ABS(aic(i)-aic(i-1))
       ENDDO
       aicder(1) = aicder(2)
       
       END
C END FILE AICD.F