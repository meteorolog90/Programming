!                                       ++++++++++++++++++++++++
!                                       ++++                ++++
!                                       ++++ ETO CALCULATOR ++++
!                                       ++++                ++++
!                                       ++++++++++++++++++++++++


!                                                                                                
!
!       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       +                                                                                 +
!       +                   0.408*delta*(Rn-G)+gama*(900/T+273)*u*(es-ea)                 +
!       +  Eto = ----------------------------------------------------------------------   +
!       +                             delta+gama*(1+0.34*u)                               +
!       +                                                                                 +
!       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!           Eto   - evapotranspiration [mm/day]
!           Rn    - global net radiation [MJ/m**2day
!           delta - slope vapour pressure [kPa]
!           G     - soil heat flux [MJ/m**2day]
!           T     - temperature at 2m height [C]
!           u     - wins speed at 2m height [m/s]
!           es    - saturation vapour pressure [kPa]
!           ea    - actual vapour pressure [kPa]
!           es-ea - saturation vapour pressure deficit [kPa]
!           gama  - psychrometric constant [kPa/C]
!      
!                       4098*(0.6108*exp(17.27*T/T+237.3)
!           delta = -----------------------------------------
!                             (T+237.3)**2
!       
!                        
!  


      program eto
      
      implicit none

      integer :: i
      real,parameter :: p=101.3,gama=p*0.66E-3
      integer,parameter :: N=9861,k=0
      integer,dimension(N):: godina,mesec,dan 
      real,dimension(N)::tmax,tmin,tsr,esr,fsr,ss,nsr,rr,delta,EKS,G,es&
      ,D,Usr,Rn,Tra,Tra1,Tra2,EtoR,DEF,es_tmax,es_tmin,esr1,tdew,ea
      
      open(1,file="Vojvodina_RHMZ_f.csv",status="old")
      
      do i=1,N
      read (1,*)godina(i),mesec(i),dan(i),tmax(i),tmin(i),tsr(i),esr(i)&
                ,fsr(i),ss(i),nsr(i),rr(i)
      
!     tmx  - maksimalna temperatrua [C]
!     tmin - minimalna temperatura [C]
!     tsr  - srednja temperatura [C]
!     esr  - srednji napon vodene pare [mb]
!     fsr  - srednja jacina vetra [bofor]
!     ss   - trajanje sijanje sunca [h]
!     nsr  - srednja oblacnost [desetine pokrivenosti neba]
!     rr   - kolicina padavina [mm]

      if (nsr(i).eq.-9999.) then
      print *, "nema podataka"
      endif
     
      if (mesec(i).EQ.1.)  Rn(i)=(920.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.2.)  Rn(i)=(1440.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.3.)  Rn(i)=(2102.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.4.)  Rn(i)=(2856.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.5.)  Rn(i)=(3383.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.6.)  Rn(i)=(3651.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.7.)  Rn(i)=(3551.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.8.)  Rn(i)=(3115.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.9.)  Rn(i)=(2437.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.10.) Rn(i)=(1717.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.11.) Rn(i)=(1097.*(1.-0.89*nsr(i)/10))/100
      if (mesec(i).EQ.12.) Rn(i)=(796.*(1.-0.89*nsr(i)/10))/100
      
      if (esr(i).LE.0.) esr(i)=esr(i)*(-1.)

      if (tsr(i).LE.0.) then 
      tdew(i)=(237.3*(alog10(esr(i))-0.786))/(7.5+0.786-alog10(esr(i)))
      else
      tdew(i)=(237.3*(alog10(esr(i))-0.786))/(9.5+0.786-alog10(esr(i)))
      endif

!     tdew - tacka rose        


      ea(i)=(0.6108*exp((17.27*tdew(i))/(tdew(i)+273.3)))/10    
      es_tmin(i)=0.6108*exp((17.27*tmin(i))/(tmin(i)+237.3))
      es_tmax(i)=0.6108*exp((17.27*tmax(i))/(tmax(i)+237.3))
      es(i)=(es_tmin(i)+es_tmax(i))/2
      EKS(i)=exp(17.27*tsr(i)/(tsr(i)+267.3))
      delta(i)=(4098.*0.6108*EKS(i))/(tsr(i)+237.3)**2 
      G(i)= Rn(i)*k
      Usr(i)=sqrt(fsr(i)**3)*0.836
      DEF(i)=es(i)-ea(i)

      Tra(i)=0.408*delta(i)*((Rn(i))-G(i))
      Tra1(i)=gama*(900./((tsr(i)+273)*Usr(i)*DEF(i)))
      Tra2(i)=delta(i)+gama*(1.+(0.34*Usr(i)))
      EtoR(i)=(Tra(i)+Tra1(i))/Tra2(i)

10    FORMAT(f6.2)
      write (*,10) EtoR(i)
      print *


      enddo
      end program eto 
