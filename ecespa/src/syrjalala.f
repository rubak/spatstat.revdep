c-------------------------------------------------------------------
      subroutine syrjalala(x,y,var1,var2,nd,nperm,cvm,ks)
      parameter (epsilon = 1.0e-16)
c
c     A program to examine a proposed statistical test for a
c     difference between two underlying geographic (or other spatial)
c     distributions.  Examined are bivariate generalizations of the
c     Cramer-von Mises test and the Kolmogorov-Smirnov test.
c     The level of significance of the test statistics is determined
c     using a randomization and/or an approximate randomization test.
c
c     The development of the statistical test is presented in the
c     paper "A Statistical Test for a Difference between the Spatial
c     Distributions of Two Populations" published in 1996, Ecology,
c     Volume 77(1):75-80.
c
c     Written by Stephen E. Syrjala                  February 1994
c     Alaska Fisheries Science Center, RACE Division
c     National Marine Fisheries Service, NOAA
c     7600 Sand Point Way Northeast
c     Bin C15700
c     Seattle, Washington 98115-0070  U.S.A.
c     e-mail:  SyrjalaS@AFSC.NOAA.GOV
c     Phone:  206-526-4135
c
c     The program on diskette is distributed as ESPS Document C9503
c     by         The Ecological Society of America
c                328 East State Street
c                Ithaca, NY 14850-4318  U.S.A.
c
c     on SAS-IML:
c     Jose M. Blanco-Moreno                     March 2006
c     Department of Biological Sciences
c     University of Alberta
c     Edmonton
c
c     on Fortran:
c     Jose M. Blanco-Moreno                     October 2008
c     Department of Plant Biology
c     University of Barcelona
c
c     it only accepts a partial randomization test; this is not a 
c     big issue, since complete randomization is feasible only for 
c     very few observations.
c
c     negative values checked for in the r function
c     duplicated coordinates are checked in the r function
c     
c     History:
c     2008-10-27  The random number generator of the permuteobs 
c     subroutine has been changed to the R uniform random number
c     generator.
c     2020-02-10 Erased unused variable 'timearray' . 
c     2020-02-10 Defined double precision variable 'halfvalue" 
c     to avoid conflicts  in comparisons with 'u' (previously 'u' 
c     was compared to 0.5, that implicitly was defined as 'real')
c
      integer  nd,nperm
      double precision  x(nd),y(nd),var1(nd),var2(nd),
     +         cvm(nperm+1),ks(nperm+1)
c      
      integer  i
      double precision  svar1,svar2,pvar1(nd),pvar2(nd),
     +         minx,miny,maxx,maxy,rangex,rangey,
     +         ppvar1(nd),ppvar2(nd),
     +         tpcvm,tpks
c
      svar1 = 0.0
      svar2 = 0.0
      minx = 1.0e25
      miny = 1.0e25
      maxx = -1.0e25
      maxy = -1.0e25
      do i = 1,nd
         svar1 = svar1 + var1(i)
         svar2 = svar2 + var2(i)
         if (x(i).lt.minx) minx = x(i)
         if (y(i).lt.miny) miny = y(i)
         if (x(i).gt.maxx) maxx = x(i)
         if (y(i).gt.maxy) maxy = y(i)
      end do
      rangex = maxx - minx
      rangey = maxy - miny
      if (rangex.le.epsilon) then rangex = dble(1.0)
      if (rangey.le.epsilon) then rangey = dble(1.0)
      do i = 1,nd
         pvar1(i) = var1(i)/dble(svar1)
         pvar2(i) = var2(i)/dble(svar2)
         x(i) = (x(i) - minx + epsilon)/rangex
         y(i) = (y(i) - miny + epsilon)/rangey
      end do
c
      call teststat(x,y,nd,pvar1,pvar2,tpcvm,tpks)
      cvm(1) = tpcvm
      ks(1) = tpks
c
c     permute and recompute the test statistic
c
      do i = 1,nperm
      	call permuteobs(pvar1,pvar2,nd,ppvar1,ppvar2)
      	call teststat(x,y,nd,ppvar1,ppvar2,tpcvm,tpks)
         cvm(i+1) = tpcvm         
         ks(i+1) = tpks
      end do
c      
      return
      end
c-------------------------------------------------------------------
c-------------------------------------------------------------------
      subroutine dcumulative(x,y,pvar1,pvar2,nd,ssdiff,maxdiff)
c
c     Subroutine to compute the cumulative distribution function of
c     data. The cumulative distribution function is defined as the
c     sum of the contribution of each datum[j] that has 
c     ((x[j]<=x[i])&(y[j]<=y[i])) (That means: is enclosed in the
c     rectangle from the origin of coordinates)
c
      integer  nd
      double precision  x(nd),y(nd),pvar1(nd),pvar2(nd),ssdiff,
     +         maxdiff
      
      integer  i,j
      double precision cumvar1,cumvar2,tmp0
c      
      maxdiff = 0.0
      ssdiff = 0.0
c      
      do i = 1,nd
         cumvar1 = 0.0
         cumvar2 = 0.0
         do j = 1,nd
            if (x(j).le.x(i).and.y(j).le.y(i)) then
               cumvar1 = cumvar1 + pvar1(j)
               cumvar2 = cumvar2 + pvar2(j)
            end if
         end do
         tmp0 = abs(cumvar1 - cumvar2)
         ssdiff = ssdiff + tmp0*tmp0
         if (tmp0.gt.maxdiff) then
            maxdiff = tmp0
         end if
      end do
c
      return
c
      end
c-------------------------------------------------------------------
c-------------------------------------------------------------------
      subroutine teststat(x,y,nd,pvar1,pvar2,cvm,ks)
c     compute the statistics for each permutation of data
      parameter(epsilon = 1.0e-16)
      integer  nd
      double precision  x(nd),y(nd),pvar1(nd),pvar2(nd),cvm,ks
c     
      integer  i      
      double precision  xx(nd),yy(nd),tssdiff,tmaxdiff
c     
      cvm = 0.0
      ks = 0.0
      tssdiff = 0.0
      tmaxdiff = 0.0
c     
c     calculate statistics from "lower left" corner
c
      do i = 1,nd
         xx(i) = x(i)
         yy(i) = y(i)
      end do
      call dcumulative(xx,yy,pvar1,pvar2,nd,tssdiff,tmaxdiff)
      cvm = cvm + tssdiff/dble(4.0)
      ks = ks + tmaxdiff/dble(4.0)
c
c     calculate statistics from "upper left"
c
      do i = 1,nd
         yy(i) = -1.0*yy(i)
      end do
      call dcumulative (xx,yy,pvar1,pvar2,nd,tssdiff,tmaxdiff);
      cvm = cvm + tssdiff/dble(4.0)
      ks = ks + tmaxdiff/dble(4.0)
c
c     calculate statistics from "upper right"
c
      do i = 1,nd
         xx(i) = -1.0*xx(i)
      end do
      call dcumulative (xx,yy,pvar1,pvar2,nd,tssdiff,tmaxdiff);
      cvm = cvm + tssdiff/dble(4.0)
      ks = ks + tmaxdiff/dble(4.0)
c
c     calculate statistics from "lower right"
c
      do i = 1,nd
         yy(i) = -1.0*yy(i)
      end do
      call dcumulative (xx,yy,pvar1,pvar2,nd,tssdiff,tmaxdiff);
      cvm = cvm + tssdiff/dble(4.0)
      ks = ks + tmaxdiff/dble(4.0)
c      
      return
      end
c-------------------------------------------------------------------
c-------------------------------------------------------------------
      subroutine permuteobs(pvar1,pvar2,nd,ppvar1,ppvar2)         
c                                                                 
c     interchange of observations among variables:                
c     this is not a "position" permutational test; what is tested 
c     instead is the difference between sampled populations, that 
c     is, if their values are exchangeable.                       
c
c     2008-10-27 This version of permuteobs is slightly slower 
c     than the original one, but it uses a call to R uniform random
c     number generator through a C wrapper function (in wrp.c)
      
      integer  nd                                                 
      double precision  pvar1(nd),pvar2(nd),ppvar1(nd),ppvar2(nd) 
                                                                  
      integer i                                   
c      double precision rand,u,sumppvar1                           
c      double precision u,sumppvar1                           
      double precision sumppvar1
      REAL u
                                                                  
      sumppvar1 = 0.0  
c                                                                 
c     inter-change values between variables if rand > 0.5         
c                                                                 
      call rndstart()
      do i = 1,nd 
          u = unifrnd()
c         u = rand(0)
         if (u.gt.0.5) then
            ppvar1(i) = pvar2(i)                                  
            ppvar2(i) = pvar1(i)                                  
         else                                                     
            ppvar1(i) = pvar1(i)                                  
            ppvar2(i) = pvar2(i)                                  
         end if                                                   
         sumppvar1 = sumppvar1 + ppvar1(i)                        
      end do                                                      
      call rndend()
c                                                                 
c     re-standardize the data in each array                       
c                                                                 
      do i = 1,nd                                                 
         ppvar1(i) = ppvar1(i)/sumppvar1                          
         ppvar2(i) = ppvar2(i)/(2.0 - sumppvar1)                  
      end do                                                      
c                                                                 
      return                                                      
      end
c-------------------------------------------------------------------
c-------------------------------------------------------------------
