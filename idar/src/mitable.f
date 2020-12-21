c------------------------------------------------------------------- 
      SUBROUTINE mtb(x1,y1,x2,y2,l1,l2,nombresp,nsp,r,nr,ltab,
     +                    abu)
      integer l1,l2,nombresp(l2),nsp,nr,ltab
      double precision x1(l1),y1(l1),x2(l2),y2(l2), r(nr),
     +                 abu(ltab*nr)

      integer i,j,nabu, rabu !,n
      double precision zero,dist  !(l1*l2)
      zero=0.0
      DO i= 1,l1
         DO j= 1,l2
c            n=n+1
            nabu =  nsp*(i-1) + nombresp(j)    ! identifica la especie en la que anotar la cuenta       
c            dist(n) = sqrt((((x1(i)-x2(j))**2)+((y1(i)-y2(j))**2)))
            dist = sqrt((((x1(i)-x2(j))**2)+((y1(i)-y2(j))**2)))
            DO k=1,nr
c               IF(dist(n).le.r(k)) THEN
               IF(dist.gt.zero.and.dist.le.r(k)) THEN
                 rabu= (l1*nsp*(k-1))+nabu   ! para la especie identificada, busca la distancia en la que anotar la cuenta
                 abu(rabu) = abu(rabu) + 1
               END IF
            END DO   
         END DO
         
         
      END DO
      
      RETURN
      END

c-------------------------------------------------------------------       
c     mtb; renombrado el nombre de la subrutina para que no conflicte 
c          con el nombre de la función en R
c     mitable: renombrada sin numero para incluir en el paquete ISAR
c       es idéntca a mitable4
c     mitable4: se diferencia de mitable3 en que incorpora en el IF
c               que "dist.gt.0", para no contar el propio arbol focal
c     mitable3: se diferencia de mitable2 en que "dist" es interna
c      al programa fortran y no se devuelve a R, con lo que en teoría
c      se debería aumantar rapidez y limmitar problemas de memoria
c           
c     mitable: calcula vecinos de cada especie  del patron 2 
c     [coordenadas (x2, y2), con l2 individuos],dentro del círculo
c     de radio r  alrededor de cada individuo del patron 1
c     [coordenadas(x1, y1), con l1 individuos]
c     Otras variables suministradas:
c     nsp: numero de especies diferentes en el patron 2.Es un integer
c            
c     nombresp: vector de integers con los numeros que representan
c     a cada especie del patron2.se obtiene con
c      unique(as.numeric(factor(marks)))
c     dist: vector de distancias que se quiere calcular
c     abu: vector de abundancias de cada especie como vecinas de 
c     cada punto del patron 1. Es la representación vectorizada de
c     la tabla de "comunidades".
 
