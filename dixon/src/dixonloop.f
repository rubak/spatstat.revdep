c------------------------------------------------------------------- 
      SUBROUTINE dixonloop(k,N,R,Q,SP,VN,VarN,EN)
c
c  Fortran code included to speed the computation of Dixon (2002) 
c  NN contingence table indices. 09-14-2010
c  There remains annotations of the original R code to facilitate debuging
c  Author: Marcelino de la Cruz, based in the R code of P.M.Dixon
c
c
      integer k, VN(k*k*k*k)
      double precision N, R, Q, SP(k), VarN(k*k*k*k), EN(k*k)

      INTEGER  l1, l2, i, j, i2, j2, kk, posij, posi2j2, posl,posl2 
      DOUBLE PRECISION p2, p3, p4, vp 
    
      kk = k*k ! para los loops
      
      DO l1= 1,kk !empieza el primer loop    
c     for (l1 in 1:(k * k)) {
        i = 1 + (l1 - 1)/k
        j = 1 + MOD(l1 - 1, k)
        DO l2= l1,kk !empieza el segundo loop    
c       for (l2 in l1:(k * k)) {
            i2 = 1 + (l2 - 1)/k
            j2 = 1 + MOD(l2 - 1, k)
            posl = ((l2-1)*k*k)+l1 !posición en las matrices vectorizadas grandes (VN[l1,l2] y VarN[l1,l2])
            posl2 = ((l1-1)*k*k)+l2 !posición en las matrices vectorizadas grandes (VN[l2,l1] y VarN[l2,l1])  
            posij =  ((j-1)*k)+i     !posición en las matrices vectorizadas pequeñas (EN[i,j])
            posi2j2 =  ((j2-1)*k)+i2  !posición en las matrices vectorizadas pequeñas (EN[i2,j2])
c            if ((i == i2) & (j == j2)) {
          IF (i.eq.i2.and.j.eq.j2.and.i.eq.j) THEN  
c                if (i == j) {
c                IF(i.eq.j) THEN
            p2 = SP(i)*(SP(i) - 1)/(N * (N - 1))
            p3 = p2 * (SP(i) - 2)/(N - 2)
            p4 = p3 * (SP(i) - 3)/(N - 3)
c                  VN[l1,l2]<- 1 # no es necesaria la función check. si el código está bien VN[l1,l2] era 0 antes
            VN(posl)=1
      vp=(N+R)*p2+(2*N-2*R+Q)*p3+(N*(N-3)-Q+R)*p4-EN(posij)*EN(posij)
            VarN(posl)=vp

c                  }
c                else {
                ELSE IF (i.eq.i2.and.j.eq.j2.and.i.ne.j) THEN  
                  p2 = SP(i) * SP(j)/(N * (N - 1))
                  p3 = p2 * (SP(i) - 1)/(N - 2)
                  p4 = p3 * (SP(j) - 1)/(N - 3)
                  
                  VN(posl) = 2
              vp=N*p2+Q*p3+(N*(N-3)-Q+R)*p4-EN(posij)*EN(posij)            
              
              VarN(posl)=vp 
c            END IF
c                }
c            }
c            else if ((i == j) & (i == i2) & (j != j2)) {
             ELSE IF(i.eq.j.and.i.eq.i2.and.j.ne.j2) THEN
                p3 = SP(i) * (SP(i) - 1) * SP(j2)/(N * (N - 1) * 
     +                  (N - 2))
                p4 = p3 * (SP(i) - 2)/(N - 3)
                VN(posl) = 3
      vp=(N-R)*p3+(N*(N-3)-Q+R)*p4
     +                                    - EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 
c            }
c            else if ((i2 == j2) & (i == i2) & (j != j2)) {
             ELSE IF(i2.eq.j2.and.i.eq.i2.and.j.ne.j2) THEN
               p3 = SP(i2) * (SP(i2) - 1) * SP(j)/(N * (N - 1) * 
     +                  (N - 2))
               p4 = p3 * (SP(i) - 2)/(N - 3)
                
              VN(posl) =3
              vp= (N-R)*p3 + (N*(N-3)-Q+R)*p4  
     +                                - EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i2 == j2) & (j == j2) & (i != i2)) {
             ELSE IF(i2.eq.j2.and.j.eq.j2.and.i.ne.i2) THEN
                p3 = SP(j) * (SP(j) - 1) * SP(i)/(N * (N - 1) * 
     +                (N - 2))
                p4 = p3 * (SP(j) - 2)/(N - 3)
                VN(posl) = 4  
                vp= (N-R+Q)*p3 + (N*(N-3)-Q+R)*p4
     +                                    - EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i == j) & (i == j2) & (i != i2)) {
             ELSE IF(i.eq.j.and.i.eq.j2.and.i.ne.i2) THEN
                p3 = SP(i) * (SP(i) - 1) * SP(i2)/(N * (N - 1) * 
     +            (N - 2))
                p4 = p3 * (SP(i) - 2)/(N - 3)
                VN(posl) = 4
c               VN(posl) = 14 !!! OJO, This looks as a typo in the original R code.Probably it meaned VN[l1,l2] <- 4 instead of <- 14
                vp= (N-R+Q)*p3 + (N *(N-3)-Q+R)*p4  
     +                                   -EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i == j) & (i2 == j2) & (i != i2)) {
             ELSE IF(i.eq.j.and.i2.eq.j2.and.i.ne.i2) THEN  
c                 
            p4 =SP(i)*(SP(i)-1)*SP(i2)*(SP(i2)-1)/(N*(N-1)*(N-2)*(N-3))
                VN(posl)= 5

                vp=  (N*(N-3)-Q+R)*p4-EN(posij)*EN(posi2j2)
             
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i == j) & (i2 != i) & (j2 != j) & (i2 != j2)) {
             ELSE IF(i.eq.j.and.i2.ne.i.and.j2.ne.j.and.i2.ne.j2) THEN 
c                              
                p4=SP(i)*(SP(i)-1)*SP(i2)*SP(j2)/(N*(N-1)*(N-2)*(N-3))
                VN(posl)= 6
                vp=  (N*(N-3)-Q+R)*p4-EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i2 == j2) & (i2 != i) & (j2 != j) & (i !=  j)) {
             ELSE IF(i2.eq.j2.and.i2.ne.i.and.j2.ne.j.and.i.ne.j) THEN  
                p4=SP(i2)*(SP(i2)-1)*SP(i)*SP(j)/(N*(N-1)*(N-2)*(N-3))
                VN(posl)= 6
                vp=  (N*(N-3)-Q+R)*p4 - EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i == i2) & (i != j) & (i2 != j2) & (j != j2)) {
             ELSE IF(i.eq.i2.and.i.ne.j.and.i2.ne.j2.and.j.ne.j2) THEN  
c                              
                p4 = SP(i) * (SP(i) - 1) * SP(j) * SP(j2)/(N * (N - 
     +            1) * (N - 2) * (N - 3))
                VN(posl)= 7  
                vp= (N*(N-3)-Q+R)*p4 - EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i == j2) & (i2 == j) & (i != j)) {
             ELSE IF(i.eq.j2.and.i2.eq.j.and.i.ne.j) THEN  
                p2 = SP(i) * SP(j)/(N * (N - 1))
                p3 = p2 * (SP(i) - 1 + SP(j) - 1)/(N - 2)
                p4 = p2 * (SP(i) - 1) * (SP(j) - 1)/((N - 2) * (N - 3))
                VN(posl)= 8             
                vp= R*p2 + (N-R)*p3 + (N*(N-3)-Q+R)*p4
     +                                   -EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i != j) & (j == i2) & (i2 != j2) & (i != j2)) {
             ELSE IF(i.ne.j.and.j.eq.i2.and.i2.ne.j2.and.i.ne.j2) THEN  
                p3 = SP(i) * SP(j) * SP(j2)/(N * (N - 1) * (N - 2))
                p4 = p3 * (SP(j) - 1)/(N - 3)
                VN(posl)= 9  
                vp= (N-R)*p3 + (N*(N-3)-Q+R)*p4
     +                                  -EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i != j) & (j == j2) & (i2 != j2) & (i != i2)) {
             ELSE IF(i.ne.j.and.j.eq.j2.and.i2.ne.j2.and.i.ne.i2) THEN  
                p3 = SP(i) * SP(j) * SP(i2)/(N * (N - 1) * (N - 2))
                p4 = p3 * (SP(j) - 1)/(N - 3)
                VN(posl)= 10  
                vp= Q*p3 + (N*(N - 3)-Q+R)*p4
     +                                   - EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i != j) & (i == j2) & (i2 != j2) & (j !=  i2)) {
             ELSE IF(i.ne.j.and.i.eq.j2.and.i2.ne.j2.and.j.ne.i2) THEN  
                p3 = SP(i) * SP(j) * SP(i2)/(N * (N - 1) * (N - 2))
                p4 = p3 * (SP(i) - 1)/(N - 3)
                VN(posl)= 11  
                vp=(N-R)*p3 +(N *(N-3)-Q+R)*p4
     +                                  -EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
c            else if ((i != j) & (i != i2) & (i != j2) & (j != i2) & (j != j2) & (i2 != j2)) {
             ELSE IF(i.ne.j.and.i.ne.i2.and.i.ne.j2.and.j.ne.i2.
     +               and.j.ne.j2.and.i2.ne.j2) THEN  
c                      
                p4= SP(i)*SP(j)*SP(i2)*SP(j2) / (N*(N-1)*(N-2)*(N-3))
                VN(posl)= 12  
                vp=(N*(N-3)-Q+R)*p4-EN(posij)*EN(posi2j2)
             VarN(posl)= vp   
             VarN(posl2)=vp 

c            }
             END IF
        END DO  !acaba el segundo loop
      END DO !acaba el primer loop 
      RETURN
      END
      
