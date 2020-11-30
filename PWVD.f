*************************************************************************
* Programa de cálculo para la PWVD                                      *
*                                                                       *
* Con opción de aplicar un filtro digital pasa-altas                    *
* Basado en el código de Shin y Jeon (1993)                             *
*************************************************************************

*************************************************************************
* Se pide el nombre de entrada de una serie de tiempo con extensión .dat*
* que consta de dos columnas. En la primera es la historia de tiempos y *
* en la segunda es la historia de las amplitudes.                       *
*                                                                       *
* Variables del programa:                                               * 
*                                                                       *
*  dp: número de datos muestreados,                                     *
*  mm: parámetro de salida de tiempo,                                   *
*  nn: parámetro de salida de frecuencia,                               *
*  mt: parámetro de la ventana Gaussiana en el tiempo,                  *
*  nf: parámetro de la ventana Gaussiana en la frecuencia,              *
*  inname: nombre de la serie de tiempo de entrada,                     *
*  mvopt: opción la cual decide si remueve o no el valor medio de la    *
*         señal,                                                        *
*  df: resolución en la frecuencia,                                     *
*  dt: intervalo de muestro.                                            *
*                                                                       *
* Arreglos de datos:                                                    *
*                                                                       *
*  tin(i): datos de tiempo muestreados,                                 *
*  ain(i): datos de la magnitud de la amplitud muestreada,              *
*  fain(i): datos filtrados,                                            *
*  wdf(i,j): arreglo reducido de la PWVD,                               *
*  bk(i): pesos del filtro,                                             *
*  s(i): señal analítica,                                               *
*  c(i): auto-correlación o datos resultantes de la FFT                 *
*************************************************************************

*************************************************************************
*                     Inicialización del programa                       *
*                                                                       *
      program PWVD
*                                                                       *
*************************************************************************
*************************************************************************
*                     Declaración de variables                          *
*                                                                       *
* -Para series de tiempo más grandes que 2048 muestras cambie np.       *
*                                                                       *
      integer dp,dp2,mvopt,redopt,mm,nn,nf,mt
      parameter (np = 4096)
      real pi,tin(np),ain(np),wdf(256,128),dt,fain(np),bk(1100)
      complex s(np*2),c(np*2)
      character*25 inname
      character as
      pi=atan(1.)*4.
*                                                                       *
*************************************************************************
*************************************************************************
*                      Impresiones de pantalla                          *
*                                                                       *
      print*,'***********************************************'
      print*,'*             Cálculo de la PWVD              *'
      print*,'***********************************************'
*                                                                       *
      print*
      print*,'***********************************************'
      print*,'Ingrese el nombre de la señal de entrada       '
      print*
      read(5,901) inname
      print*
      print*,'***********************************************'
      print*,'¿Cuál es el número total de muestras?          '
      print*
      read(5,902) dp
      print*
*                                                                       *
      dp2=dp*2
*                                                                       *
      do 100 i=1,dp2
 100     s(i)=cmplx(0.,0.)
*                                                                       *
      print*,'***********************************************'
      print*,'Desea remover el valor medio de la señal       '
      print*
      print*,'1 para sí, 0 para no                           '
      print*
      read(5,902) mvopt
      print*
*                                                                       *
      print*,'***********************************************'
      print*,'¿Desea aplicar el filtro digital? (Y/N)        '
      print*
      read(*,903) as
      print*
*                                                                       *
      if(as.eq.'Y'.or.as.eq.'y') then
         print*,'***********************************************'
         print*,'Ingrese la frecuencia de corte (Hz)            '
         print*
         read(*,*) bw
         print*
      endif
*                                                                       *
      fmin = 0.
*                                                                       *
      print*,'***********************************************'
      print*,'Ingrese el tamaño de salida de la PWVD en      '
      print*,'Tiempo vs. Frecuencia                          '
      print*,'_______________________________________________'
      print*,'Ingrese 1 para 64 por 32                       '
      print*,'Ingrese 2 para 128 por 64                      '
      print*,'Ingrese 3 para 128 por 128                     '
      print*,'Ingrese 4 para 256 por 128                     '
      print*
      read(5,902) redopt
      print*
      if(redopt.eq.1) then
         mm = dp/32
         nn = dp2/64
      elseif(redopt.eq.2) then
         mm = dp/64
         nn = dp2/128
      elseif(redopt.eq.3) then
         mm = dp/128
         nn = dp2/128
      else
         mm = dp/128
         nn = dp2/256
      endif
*                                                                       *
      print*,'***********************************************'
      print*,'Ingrese los parámetros de la ventana Gaussiana '
      print*
      print*,'Parámetro de frecuencia:                       '
      print*
      read(5,902) nf
      print*
      print*,'Parámetro de tiempo:                           '
      print*
      read(5,902) mt
      print*
*                                                                       *
*************************************************************************
*************************************************************************
*                      Cálculos                                         *
*                                                                       *
* Leé los datos de la señal de entrada:
*
      open(4, file = inname, status = 'old')
      call indata(dp, tin, ain)
      print*,'***********************************************'
      print*,'* Lectura de datos completada                  '
      print*
      close(4)
*                                                                       *
* Calcula el incremento del tiempo de muestreo:                         *
*                                                                       *
      call dtcalc(dp, tin, dt)
      print*,'* Cálculo de la delta de tiempo finalizado     '
      print*
*                                                                       *
* Calcula el valor medio de la señal:                                   *
*                                                                       *
      if(mvopt.eq.1) then
         call mean(dp, ain)
         print*,'* Cálculo del valor medio terminado         '
         print*
      endif
*                                                                       *
* -Modificaciones a la señal-                                           *
*                                                                       *
* Aplicación del filtro digital:                                        *
*                                                                       *
      if(as.eq.'Y'.or.as.eq.'y') then
         mo = dp/2
*                                                                       *
         call filter(mo, bw, dt, bk)
*                                                                       *
         do 160 i = 1,dp
 160        fain(i) = 0.
*                                                                       *
         do 200 i = 1,dp
            do 170 k = -mo,mo
               j=k
               if(k.lt.0) then
                  j=k*(-1)
               endif
               j = j+1
               ll = i-k
               if(ll.lt.1) then
                  ll = ll+dp
               elseif(ll.gt.dp) then
                  ll = ll-dp
               endif
               bb = -bk(j)
               if(k.eq.0) then
                  bb = 1-bk(j)
               endif
               fain(i) = fain(i)+bb*ain(ll)
 170        continue
 200     continue
*                                                                       *
         do 210 i = 1,dp
 210        ain(i) = fain(i)
*                                                                       *
        print*,'* Aplicación del filtro digital terminado   '
        print*
        endif
*                                                                       *
*                                                                       *
* Aplicación de la ventana de Hamming                                   *
*                                                                       *
      call hammg(dp, dt, pi, ain)
      print*,'* Aplicación de la ventana de Hamming terminada'
      print*
*                                                                       *
* Conversión a una señal analítca                                       *
*                                                                       *
      call anal(dp, pi, ain, s)
      print*,'* Terminada la conversión a señal analítica    '
      print*
* 
* Escribe los archivos de salida de las distribuciones                  *
*                                                                       *
      open(9, file = 'rwdf.out', status = 'new')
      open(10, file  = 'rswdf.out', status = 'new')
*
* Calcúla los parámetros de salida                                      *
*                                                                       *
      ttime = dp*dt
      df = 1./(4.*dp*dt)
      fmax = 2.*dp*df
      nx = dp2/nn
      ny = dp/mm
      write(9,904) tin(1), ttime, fmin, fmax
      write(9,*) nx, ny
*                                                                       *
* Cálculo de la PWVD                                                    *
*                                                                       *
      print*,'***********************************************'
      print*,'*              Calculando PWVD...             *'
      print*,'***********************************************'
      print*
      call wigner(dp, dt, pi, s, c, mm, nn, wdf, mt, nf)
*                                                                       *
      print*,'***********************************************'
      print*,'*              Cálculo finalizado             *'
      print*,'***********************************************'
      print*
      close(9)
*                                                                       *
* Escribe los resutados                                                 *
*                                                                       *
      write(10,904) tin(1), ttime, fmin, fmax
      write(10,*) nx, ny
      do 500 i = 1,dp/mm
         do 500 j = 1,dp2/nn
            write(10,905) wdf(j,i)
 500     continue
*                                                                       *     
901   format(a25)
902   format(i6)
903   format(a1)
904   format(e12.5,/,e12.5,/,e12.5,/,e12.5)
905   format(2x,e12.5)
*                                                                       *
      close(10)
      print*,'***********************************************'
      print*,'*      ¡Éxito! Ejecución de PVWD terminada    *'
      print*,'***********************************************'
      print*
*                                                                       *
      stop
      end
*                                                                       *
*************************************************************************
*                   Subrutinas de ejecución                             *
*************************************************************************
*                                                                       *
* Subrutina de ingreso de datos                                         *
*                                                                       *
      subroutine indata(dp, tin, ain)
*                                                                       *
      integer dp
      real tin(*), ain(*)
*                                                                       *
      do 100 j = 1,dp
         read(4,*) tin(j), ain(j)
 100  continue
      return
      end
*                                                                       *
* Subrutina de cálculo del incremento de tiempo
*                                                                       *
      subroutine dtcalc(dp, tin, dt)
*                                                                       *
      integer dp
      real tin(*), dt
*                                                                       *
      dtsum = 0.0
*                                                                       *
      do 100 i = 1,dp-1
         delt = tin(i+1) - tin(i)
         dtsum = dtsum + delt
 100  continue
*                                                                       *
      dt = dtsum/float(dp-1)
*                                                                       *
      return
      end
*                                                                       *
* Subrutina de cálculo de la media de los datos                         *
*                                                                       *
      subroutine mean(dp, ain)
*                                                                       *
      integer dp
      real ain(*), meanv
*                                                                       *
      asum = 0.0
*                                                                       *
      do 100 i = 1,dp
         asum = asum + ain(i)
 100  continue
*                                                                       *
      meanv = asum/dp
*                                                                       *
      do 200 i = 1,dp
         ain(i) = ain(i) - meanv
 200  continue
*                                                                       *
      return
      end
*                                                                       *
* Subrutina de cálculo de elementos del filtro digital                  *
*                                                                       *
      subroutine filter(mo, b, t, bk)
*                                                                       *
      dimension bk(*),d(3)
      data d0/0.35577019/,d(1)/0.2436983/,d(2)/0.07211497/, 
     *d(3)/0.00630165/
*                                                                       *
      pi = atan(1.)*4
      m = mo
      fact = 2.*b*t
      bk(1) = fact
      fact = fact*pi
*                                                                       *      
      do 5 i = 1,m
         fi = i
 5       bk(i+1) = sin(fact*fi)/(pi*fi)
*                                                                       *
      bk(m+1) = bk(m+1)/2.
      sumg = bk(1)
*                                                                       *      
      do 15 i = 1,m
         sum = d0
         fact = pi*float(i)/float(m)
         do 10 k = 1,3
 10         sum = sum + 2.*d(k)*cos(fact*float(k))
         bk(i+1) = bk(i+1)*sum
 15      sumg = sumg + 2.*bk(i+1)
*                                                                       * 
         m1 = m+1
*                                                                       *         
         do 20 i = 1,m1
 20         bk(i) = bk(i)/sumg
*                                                                       *
         return
         end
*                                                                       *
* Subrutina del cálculo de la ventana de Hamming                        *
*                                                                       *
      subroutine hammg(dp, dt, pi, ain)
*                                                                       *
      integer dp
      real pi, ain(*), dt, mtime, del1, del2, const
*                                                                       *
      mtime = (dp-1)*dt
      del1 = 0.1*mtime
      del2 = 0.9*mtime
      const = pi/del1
*                                                                       *
      do 100 j = 1,dp
         t = (j-1)*dt
         if(t.le.del1) then
            ain(j) = ain(j)*(.54-.46*cos(const*t))
         elseif((t.ge.del2).and.(t.le.mtime)) then
            ain(j) = ain(j)*(.54-.46*cos(const*(mtime-t)))
         endif
 100  continue
*                                                                       *      
      return
      end
*                                                                       *
* Subrutina de conversión a señal analítica
*                                                                       *
      subroutine anal(dp, pi, ain, s)
*                                                                       *
      integer dp
      real pi, ain(*), sum, sumb, val, sval
      complex s(*)
*                                                                       *
      do 100 i = 1,dp
         sum = 0.0
         do 200 j = 1,dp
            sumb = 0.0
            if(i-j.eq.0) go to 200                                            
            n = i-j
            val = pi*n/2.
            sval = sin(val)
            sumb = ain(j)*sval*sval/val
 200     sum = sum+sumb
      s(i) = cmplx(ain(i), sum)
 100  continue
*                                                                        *
      return
      end
*                                                                        *
* Subrutina de cálculo de la distribución de Wigner
*                                                                        *
      subroutine wigner(dp, dt, pi, s, c, mm, nn, wdf, mt, nf)
*                                                                        *
      integer dp, dp2
      real pi, dt, coef, wdf(256,128), hg(-100:100,-100:100)
      complex s(*), dum, c(*)
*                                                                        *
      dp2 = dp*2
      coef = 2.0*dt
      df = 1./(4.*dp*dt)
      nf2 = nf*2
      mt2 = mt*2
      f1 = float(mt)
      f2 = float(nf)
*
*  Cálculo anidado de la función de ventana Gaussiana                     *
*      
      val = 1./((2.*pi)**2*f1*f2*df*dt)
*                                                                         *
      do 20 j = -mt2,mt2
         q1 = float(j)
         do 10 i = -nf2,nf2
            q2 = float(i)
            cf = -((q1*q1)/(2.*f1*f1))-((q2*q2)/(2.*f2*f2))
            hg(i,j) = val*exp(cf)
 10      continue
 20   continue
*                                                                         *
      do 100 i = 1,256
         do 100 j = 1,128
 100        wdf(i,j) = 0.
*                                                                         *
      do 5000 j = 1,dp
         do 1000 i = 1,dp+1
            if(j.ge.i) then
               dum = s(j-i+1)
            else
               dum = cmplx(0.,0.)
            endif
            c(i) = coef*(s(j+i-1)*conjg(dum))
            if(i.ne.1.and.i.ne.dp+1) then
               c(dp2-i+2) = conjg(c(i))
            endif
            open(11, file = 'corr.out', status = 'new')
         write(11,1400) real(c(i))
 1000    continue        
*                                                                          * 
      call fft(dp, pi, c)
      ik = mod((j-1), mm)
      if(ik.eq.0) then
         do 1500 i = 1,dp2,nn
            write(9,1400) real(c(i))
 1400       format(2x,e12.5)
 1500    continue
      endif
*                                                                           *
      m1 = 0
*                                                                           *
      do 4000 m = 1,dp,mm
         m1 = m1+1
         n1 = 0
         if(abs(j-m).le.mt2) then
            do 3000 n = 1,dp2,nn
               n1 = n1+1
               do 2500 kk = n-nf2,n+nf2
                  k1 = kk
                  if(kk.lt.1) k1 = kk+dp2
                  if(kk.gt.dp2) k1 = kk-dp2
                  wdf(n1,m1) = wdf(n1,m1)+real(c(k1))*hg(kk-n,j-m)*df*dt
 2500          continue
 3000       continue
          endif
 4000  continue
 5000  continue
       return
       end
*                                                                            *
* Subrutina de cálculo de la transformada rápida de Fourier                  *
*                                                                            *
      subroutine fft(dp, pi, c)
*                                                                            *
      integer dp, dp2, val, coef, coef1
      real pi
      complex dum, c(*), dum3, dum2
*                                                                            *
      dp2 = dp*2
      const = float(dp2)
      val = alog(const)/alog(2.)+.1
      j=1
*                                                                            *
      do 40 i = 1,dp2-1
         if(i.ge.j) go to 10
         dum3 = c(j)
         c(j) = c(i)
         c(i) = dum3
 10   k = dp
 20   if(k.ge.j) go to 30
      j = j-k
      k=k/2
      go to 20
 30   j = j+k
 40   continue
      do 70 n = 1,val
         coef = 2**n
         coef1 = coef/2
         dum2 = cmplx(1.,0.)
         theta = pi/float(coef1)
         dum = cmplx(cos(theta), -sin(theta))
         do 60 j = 1,coef1
            do 50 i = j,dp2,coef
               ii = i+coef1
               dum3 = c(ii)*dum2
               c(ii) = c(i)-dum3
               c(i) = c(i)+dum3
 50         continue
            dum2 = dum2*dum
 60      continue
 70   continue
      return
      end
*                                                                           *
*****************************************************************************
