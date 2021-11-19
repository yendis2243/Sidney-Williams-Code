!--------------------------------------------------------------------
! Программа расчёта несжимаемого ламинарного динамического
! пограничного слоя на плоской пластине по неявной схеме
! с устранением нелинейности методом простых итераций.
!--------------------------------------------------------------------
program BLIT
implicit none
!--------------------------------------------------------------------
! Блок описания переменных
real*8,   parameter      :: nu=1.8d-5   ! Кинематическая вязкость
real*8                   :: dx,dym,dyp,Kgrad,Rel,Relit,SoundSpeed,Ma ! Шаг сетки по осям x,y
integer*4,parameter      :: n=1050,k=300 ! Число узлов по осям x,y
real*8,   dimension(n,k) :: U,V,T,rho,xk,ktplus,Uplus,yplus         ! Поля скоростей u,v
integer*4                :: i,j,iglob,Turbomodel         ! Индексы текущего узла
real*8,   dimension(k,4) :: abcd  ! Трехдиагональная матрица
real*8,   dimension(n)   :: dy1,Uk,Pressure,tauw,muw,rhow  ! Skin Friction and Reynolds #
real*8                   :: M,as,b,bs,Rex0,Blas   ! Spacing constant, variable
real*8,    dimension(n,k):: y        ! Normal direction Parameter
real*8,    dimension(k)  :: dx1,lam,mu,Cp
real*8,    dimension(2,k):: x,mut
character (60)::filenameProfiles
!--------------------------------------------------------------------
! Блок итерационных переменных
real*8,   dimension(n,k) :: Uit, Vit, Tit, omegait, ktit, gamit  ! Поля итерационных
                                           ! переменных U,V
real*8                   :: epsitU, epsitV, epsitT, epsitohm, epsitkt, epsitgam ! Погрешность определения
                                           ! скоростей U,V
real*8,   parameter      :: eps = 1d-3     ! Максимальная погрешность
integer*4                :: it             ! Номер итерации
real*8                   :: maxeps         ! Функция расчёта
                                           ! максимальной погрешности
real*8,   dimension(7,k)   :: recordU,recordV,recordT,recordmu,recordmut,recordkt,recordgam,recordomega,recordRex
real*8                   :: TestRex,Rexx,Cf,Rex,St
!--------------------------------------------------------------------
! Блок расчета параметров сетки
  real*8 :: L,p,U0,del,T0,Tw,MWght,Rg,denmu,a1,a2,a3,a4,a5,a6,a7,Amu,Bmu,Cmu,Dmu,Alam,Blam,Clam,Dlam
! Turbulence Model
  real*8, dimension(n,k)              :: omega,kt,gam
  real*8                              :: ohm,pk,pgam,S,fgam,rnu,rc,egam,ggam,fturb,tomega,rt
  real*8                              :: c_mu,comega1,comega2,sigmak,sigmaomega,sigmal,sigmagamma,gammamax,Ret,Dq,TI,Tu,dyp1,dymk
!---------------------------------------------------------------------
  open (10,file='Integral Parameters.txt')
  write (10,90)'Rex','Cf','St','Blas'
  90 format(4(2x,A16))
  Turbomodel=0 !0 off, 1 on
  Rex0=1d0
  L=10d0
  p=3d0
  SoundSpeed=0d0
  Ma=0d0
  Cp(:)=0d0
  Kgrad=0d-6
  U0=5d0
  T0=300d0
  Tw=305d0
  M=1.02d0
  b=1.00003d0
  bs=0d0
  Pressure(:)=0.1d6
  MWght=28.9651784d0
  Rg=8314.15d0
  lam(:)=0.026d0
  mu(:)=2.178d-5
  mut(:,:)=0d0
  denmu=0d0
  Uk(:)=0d0
  Uk(1)=U0
  do i=0,n-2
    bs=bs+b**i
  end do
  dx1=L/bs
  del=p*L*5d0/dsqrt(U0*L/nu)
  as=0d0
  do i=0,k-2
    as=as+M**i
  end do
  dy1(:)=del/as
  y=0d0
  do j=2,k
    y(:,j)=y(:,j-1)+M**(j-2)*dy1(:)
  end do
  x(1,:)=1d-7
  xk(1,:)=0d0
  c_mu=0.09d0
  comega1=5d0/9d0
  comega2=3d0/40d0
  sigmak=2d0
  sigmaomega=2d0
  sigmal=5d0
  sigmagamma=0.2d0
  Tu=0.065d0
  TI=100d0
!--------------------------------------------------------------------
! Блок начальных данных и граничных условий
  !Initializing Arrays
  kt(:,:)=0d0
  omega(:,:)=0d0
  gam(:,:)=0d0
  U(:,:)=0d0         ! Обнуление массива для U
  U(1,:)=U0          ! U на входе
  U(:,1)=0d0         ! U на стенке
  U(1,k)=U0          ! U на внешней границе
  do j=2,k-1
    U(1,j)=U0*derf(0.313d0*y(1,j)*dsqrt(U0/(8d-4)/nu))
  end do
  V(:,:)=0d0         ! Обнуление массива для V
  V(1,:)=0d0         ! V на входе
  V(:,1)=0d0         ! V на стенке
  T(:,:)=0d0
  T(1,:)=T0
  T(:,1)=Tw
  T(:,k)=T0
  rho(:,:)=0d0
  rho(1,:)=Pressure(1)*MWght/Rg/T0
  rho(:,1)=Pressure(:)*MWght/Rg/Tw
  rho(:,k)=Pressure(:)*MWght/Rg/T0
!--------------------------------------------------------------------
!Initial Air Properties Calculation
do j=1,k
   if((50<=T(1,j)).and.(T(1,j)<350))then
    a1=-0.119634863D+03
    a2=0.700560443D+01
    a3=0.333538993D+01
    a4=0.167560348D-02
    a5=-0.906735817D-05
    a6=0.229118610D-07
    a7=-0.198110689D-10
    Amu=0.53518000E+00
    Bmu=-0.88767660E+02
    Cmu=0.24615856E+04
    Dmu=0.24432839E+01
    Alam=0.41176132E+00
    Blam=-0.15401721E+03
    Clam=0.48000858E+04
    Dlam=0.36757372E+01
   end if
   if((350<=T(1,j)).and.(T(1,j)<1000))then
    a1=1.009950160D+04
    a2=-1.968275610D+02
    a3=5.012355110D+00
    a4=-5.761013730D-03
    a5=1.066859930D-05
    a6=-7.940297970D-09
    a7=2.185231910D-12
   end if
   if((350<=T(1,j)).and.(T(1,j)<1200))then
    Amu=0.61011000E+00
    Bmu=-0.50157820E+02
    Cmu=0.10342746E+04
    Dmu=0.19056867E+01
    Alam=0.83731000E+00
    Blam=0.85130670E+02
    Clam=-0.10985937E+05
    Dlam=0.62848740E+00
   end if
   if((1000<=T(1,j)).and.(T(1,j)<5000))then
    a1=2.415214430D+05
    a2=-1.257874600D+03
    a3=5.147758670D+00
    a4=-2.138541790D-04
    a5=7.065227840D-08
    a6=-1.071483490D-11
    a7=6.577800150D-16
   end if
   if((1200<=T(1,j)).and.(T(1,j)<5000))then
    Amu=0.83884000E+00
    Bmu=0.47166816E+03
    Cmu=-0.14689912E+06
    Dmu=-0.04815000E+00
    Alam=0.88945000E+00
    Blam=0.16879525E+03
    Clam=-0.26494289E+05
    Dlam=0.19986000E+00
   end if
   Cp(j)=(a1*T(1,j)**(-2d0)+a2*T(1,j)**(-1d0)+a3+a4*T(1,j)+a5*T(1,j)**2d0+a6*T(1,j)**3d0+a7*T(1,j)**4d0)*Rg/MWght*1d+3
   mu(j)=exp(Amu*log(T(1,j))+Bmu/T(1,j)+Cmu/T(1,j)**2d0+Dmu)*1d-7
   lam(j)=exp(Alam*log(T(1,j))+Blam/T(1,j)+Clam/T(1,j)**2d0+Dlam)*1d-4
end do

if(Turbomodel==1)then
    mut(1,1)=mu(1)*TI
  kt(1,1)=1.5d0*(U(1,k)*Tu)**2d0
  gam(1,1)=0d0
  omega(1,1)=6d0*mu(1)/(comega2*rho(1,1)*(y(1,2))**(2d0))

  do j=2,k
    mut(1,j)=mu(j)*TI
    kt(1,j)=1.5d0*(U(1,k)*Tu)**2d0
    Ret=rho(1,j)/mu(j)*dsqrt(kt(1,j))*y(1,j)
    Dq=1d0-dexp(-0.022d0*Ret)
    omega(1,j)=0.09d0*Dq*rho(1,j)*kt(1,j)/mut(1,j)
    if(y(1,j)<=1e-4)then
        gam(1,j)=0d0
    else
        gam(1,j)=1d0
    end if
  end do
end if

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! Блок начальных данных для итерационных переменных
  Uit = 0d0
  Vit = 0d0
  Tit = 0d0
  omegait = 0d0
  ktit = 0d0
  gamit = 0d0
  it = 0
  iglob=1
!--------------------------------------------------------------------
! Основной цикл расчета
  do while(x(2,1)<=L)            ! Перебор по x
    i=2
    iglob=iglob+1
    rho(i,:)=rho(i-1,:)
    omega(i,:)=omega(i-1,:)
    kt(i,:)=kt(i-1,:)
    gam(i,:)=gam(i-1,:)
    x(2,:)=x(1,:)+dx1(:)*b**(iglob-1)
    xk(iglob,:)=x(2,:)
    dx=x(2,1)-x(1,1)
    dx=5d-6
    x(2,:)=x(1,:)+dx
    ! ---------------------------------------------------------------
    ! Итерационный цикл
    denmu=denmu+rho(i,k)/mu(k)*dx
    U(i,k)=-1d0/(Kgrad*denmu-1d0/U0)
    Uk(iglob)=U(i,k)
    do it =1,4000
      Uit(i,:) = U(i,:) ! Обновление итерационных переменных по U
      Vit(i,:) = V(i,:) ! Обновление итерационных переменных по V
      Tit(i,:)=T(i,:)
      omegait(i,:)=omega(i,:)
      ktit(i,:)=kt(i,:)
      gamit(i,:)=gam(i,:)
      ! ---------------------------------------------------------------
      ! Расчёт профиля скорости U
      Pressure(i)=Pressure(i-1)-rho(i,k)**2d0/mu(k)*U(i,k)**3d0*Kgrad*dx
      do j=2,k-1       ! Перебор по y
        dym=y(i,j)-y(i,j-1)
        dyp=y(i,j+1)-y(i,j)
        abcd(j,1)=-rho(i,j)*V(i,j)*dyp/dym/(dym+dyp)-2d0*(mu(j-1)+mut(i,j-1))/dym/(dym+dyp)
        abcd(j,2)= rho(i,j)*U(i,j)/dx+rho(i,j)*V(i,j)*(dyp-dym)/(dym)/dyp+2d0/(dym+dyp)*((mu(j+1)+mut(i,j+1))/dyp &
        +(mu(j-1)+mut(i,j-1))/dym)
        abcd(j,3)= rho(i,j)*V(i,j)*dym/dyp/(dym+dyp)-2d0*(mu(j+1)+mut(i,j+1))/dyp/(dym+dyp)
        abcd(j,4)= rho(i,j)*U(i,j)*U(i-1,j)/dx+rho(i-1,j)*U(i,k)*(U(i,k)-U(i-1,k))/dx
      end do
      abcd(1,1)=0d0
      abcd(1,2)=1d0
      abcd(1,3)=0d0
      abcd(1,4)=U(i,1)

      abcd(k,1)=0d0
      abcd(k,2)=1d0
      abcd(k,3)=0d0
      abcd(k,4)=U(i,k)

      call TDM(k,abcd)
      Rel=1d0
      Relit=1d0-Rel
      U(i,:)=Rel*abcd(:,4)+Relit*Uit(i,:)
       do j=2,k
        dym=y(i,j)-y(i,j-1)
        V(i,j)=rho(i,j-1)/rho(i,j)*V(i,j-1)-dym/(rho(i,j))*dx*(rho(i,j)*U(i,j)-rho(i,j)*U(i-1,j))
      end do

      do j=2,k-1
        dym=y(i,j)-y(i,j-1)
        dyp=y(i,j+1)-y(i,j)
        abcd(j,1)=-Cp(j)*rho(i,j)*V(i,j)*dyp/dym/(dym+dyp)-2d0*lam(j-1)/dym/(dym+dyp)
        abcd(j,2)=Cp(j)*rho(i,j)*U(i,j)/dx+Cp(j)*rho(i,j)*V(i,j)*(dyp-dym)/(dym)/dyp+2d0/(dym+dyp)*(lam(j+1)/dyp+lam(j-1)/dym)
        abcd(j,3)=Cp(j)*rho(i,j)*V(i,j)*dym/dyp/(dym+dyp)-2d0*lam(j+1)/dyp/(dym+dyp)
        abcd(j,4)=U(i,j)*T(i-1,j)*Cp(j)*rho(i,j)/dx
      end do
      ! На стенке j=1

      abcd(1,1)=0d0
      abcd(1,2)=1d0
      abcd(1,3)=0d0
      abcd(1,4)=T(i,1)

      ! На внешней границе j=k
      abcd(k,1)=0d0
      abcd(k,2)=1d0
      abcd(k,3)=0d0
      abcd(k,4)=T(i,k)

      call TDM(k,abcd)
      Rel=1d0
      Relit=1d0-Rel
      T(i,:)=Rel*abcd(:,4)+Relit*Tit(i,:)
if(Turbomodel==1)then
    call k_omega_gamma(dx,k,n,i,y,V,U,rho,gam,gamit,kt,ktit,omega,omegait,mu,mut)
end if


rho(i,:)=Pressure(i)*MWght/Rg/T(i,:)
      ! ---------------------------------------------------------------
      do j=1,k
        if((50<=T(i,j)).and.(T(i,j)<350))then
         a1=-0.119634863D+03
         a2=0.700560443D+01
         a3=0.333538993D+01
         a4=0.167560348D-02
         a5=-0.906735817D-05
         a6=0.229118610D-07
         a7=-0.198110689D-10
         Amu=0.53518000E+00
         Bmu=-0.88767660E+02
         Cmu=0.24615856E+04
         Dmu=0.24432839E+01
         Alam=0.41176132E+00
         Blam=-0.15401721E+03
         Clam=0.48000858E+04
         Dlam=0.36757372E+01
        end if
        if((350<=T(i,j)).and.(T(i,j)<1000))then
         a1=1.009950160D+04
         a2=-1.968275610D+02
         a3=5.012355110D+00
         a4=-5.761013730D-03
         a5=1.066859930D-05
         a6=-7.940297970D-09
         a7=2.185231910D-12
        end if
        if((350<=T(i,j)).and.(T(i,j)<1200))then
         Amu=0.61011000E+00
         Bmu=-0.50157820E+02
         Cmu=0.10342746E+04
         Dmu=0.19056867E+01
         Alam=0.83731000E+00
         Blam=0.85130670E+02
         Clam=-0.10985937E+05
         Dlam=0.62848740E+00
        end if
        if((1000<=T(i,j)).and.(T(i,j)<5000))then
         a1=2.415214430D+05
         a2=-1.257874600D+03
         a3=5.147758670D+00
         a4=-2.138541790D-04
         a5=7.065227840D-08
         a6=-1.071483490D-11
         a7=6.577800150D-16
        end if
        if((1200<=T(i,j)).and.(T(i,j)<5000))then
         Amu=0.83884000E+00
         Bmu=0.47166816E+03
         Cmu=-0.14689912E+06
         Dmu=-0.04815000E+00
         Alam=0.88945000E+00
         Blam=0.16879525E+03
         Clam=-0.26494289E+05
         Dlam=0.19986000E+00
        end if
        Cp(j)=(a1*T(i,j)**(-2d0)+a2*T(i,j)**(-1d0)+a3+a4*T(i,j)+a5*T(i,j)**2d0+a6*T(i,j)**3d0 &
        +a7*T(i,j)**4d0)*Rg/MWght
        mu(j)=exp(Amu*log(T(i,j))+Bmu/T(i,j)+Cmu/T(i,j)**2d0+Dmu)*1d-7
        lam(j)=exp(Alam*log(T(i,j))+Blam/T(i,j)+Clam/T(i,j)**2d0+Dlam)*1d-4
    end do
      ! ---------------------------------------------------------------
      ! Расчёт профиля скорости V
      ! ---------------------------------------------------------------
      epsitU = maxeps(U,Uit,n,k,i)
      epsitV = maxeps(V,Vit,n,k,i)
      epsitT = maxeps(T,Tit,n,k,i)
      epsitohm = maxeps(omega,omegait,n,k,i)
      epsitkt = maxeps(kt,ktit,n,k,i)
      epsitgam = maxeps(gam,gamit,n,k,i)
      if ((epsitU <= eps).and.(epsitV <= eps).and.(epsitT <= eps).and.(epsitohm <= eps)&
      .and.(epsitkt <= eps).and.(epsitgam <= eps)) exit
    end do
    print *, x(i,k), it
    If (it>=4000)   Print*,'epsitu= ', epsitu, 'epsitv= ', epsitv,'epsitT= ', epsitT,'epsitk= ',&
        epsitkt, 'epsitw= ', epsitohm,'epsitgam= ', epsitgam  !,epsitf, x(i,k) !,epsite,epsitk
!        If (it>=400) call Sound
!        If (it>=4000) pause
    Rex=U(i,k)*iglob*dx/(mu(k)/rho(i,k))
    Cf=mu(1)/(rho(i,1)*U(i,k)**2d0)*(U(i,2)-U(i,1))/(y(i,2)-y(i,1))
    St=(lam(1)*(T(i,2)-T(i,1))/(y(i,2)-y(i,1)))/(rho(i,1)*cp(1)*(T(i,k)-T(i,1)))
    Blas=0.332/dsqrt(Rex)
    if(mod(iglob,10)==0) write (10,100)Rex, Cf, St, mut(i,k)
    100 format(4(2x,E16.8))
    tauw(i)=(U(i,2)-U(i,1))/(y(i,2)-y(i,1))
    muw(i)=mu(1)
    rhow(i)=rho(i,1)
    do j=1,k
        yplus(iglob,j)=y(iglob,j)*dsqrt(rhow(i)*tauw(i))/muw(i)
        Uplus(i,j)=U(i,j)/dsqrt(tauw(i)/rhow(i))
        ktplus(i,j)=kt(i,j)/(tauw(i)/rhow(i))
    end do
    Rexx=rho(i,k)*U(i,k)*x(i,k)/mu(k)
    if(Rexx>=Rex0)then

        print *, Rex0
        Rex0=Rex0*10d0
        Write(filenameProfiles,'(A4,ES7.1,1X,A2,ES8.2,1X,A3,F5.0,A12)')'Rex=',Rexx,'K=',Kgrad,'Tw=',Tw,'.txt'
        open (3, file= filenameProfiles)
        write(3,70) 'y','U','V','T','rho','mu','Cp','lam'
        70 format(8(2x,A16))
        do j=1,k
            write(3,101) y(i,j),U(i,j),V(i,j),T(i,j),rho(i,j),mu(j),Cp(j),lam(j)
            101 format(8(2x,E16.8))
        end do
    end if
    U(i-1,:) = U(i,:) ! Первое приближение по U
    V(i-1,:) = V(i,:) ! Первое приближение по V
    T(i-1,:)=T(i,:)
    x(1,:)=x(2,:)
    mut(i-1,:)=mut(i,:)
    Pressure(i-1)=Pressure(i)
    SoundSpeed=331d0+0.6d0*(T(i,k)-273.15d0)
    Ma=U(i,k)/SoundSpeed
    if ((Ma>=0.4d0)) exit
    !if ((U(i,k)>=100d0)) exit
    ! ---------------------------------------------------------------
  end do
  !print *, iglob
  !pause
!--------------------------------------------------------------------
! Вывод данных о профилях скорости U,V в конце пластины
  !open (10,file='Integral Parameters.txt')
  !write (10,90)'Rex','Cf','St'
  !90 format(3(2x,A16))
  !do i=1,n
  !  write (10,100)Rex(i), Cf(i), St(i)
  !  100 format(3(2x,E16.8))
  !end do
!--------------------------------------------------------------------
pause
end program BLIT

!--------------------------------------------------------------------
! Расчёт абсолютной погрешности вычисления величины
!--------------------------------------------------------------------
function maxeps(v,vit,n,k,i)
!--------------------------------------------------------------------
! Блок описания переменных
real*8                              :: maxeps
real*8, dimension (n,k), intent(in) :: v,vit
integer*4,               intent(in) :: n,k
real*8                              :: epsi
integer*4                           :: i
!--------------------------------------------------------------------
! Основной цикл
  maxeps = dabs(v(i,1)-vit(i,1))
  do j=2,k
    epsi = dabs(v(i,j)-vit(i,j))

    if(epsi >= maxeps)then
      maxeps = epsi
    end if
  end do
!--------------------------------------------------------------------
end function maxeps

!--------------------------------------------------------------------
! Расчёт относительной погрешности вычисления величины
!--------------------------------------------------------------------
function maxepsr(v,vit,n,k,i)
!--------------------------------------------------------------------
! Блок описания переменных
real*8                              :: maxepsr
real*8, dimension (n,k), intent(in) :: v,vit
integer*4,               intent(in) :: n,k
real*8                              :: epsi
integer*4                           :: i
!--------------------------------------------------------------------
! Основной цикл
  maxepsr = dabs(v(i,1)-vit(i,1))/(vit(i,1)+1d-30)
  do j=2,k
    epsi = dabs(v(i,j)-vit(i,j))/(vit(i,j)+1d-30)
    if(epsi >= maxepsr)then
      maxepsr = epsi
    end if
  end do
!--------------------------------------------------------------------
end function maxepsr

!--------------------------------------------------------------------
! Подпрограмма расчёта трёхдиагональной матрицы
!--------------------------------------------------------------------
subroutine TDM(n,T)
!--------------------------------------------------------------------
! Блок описания переменных
integer*4, intent(in)                 :: n ! Число строк
real*8, dimension(n,4), intent(inout) :: T ! Трёхдиагональная матрица
integer*4                             :: i ! Индекс текущей строки
!--------------------------------------------------------------------
! Основные циклы
  do i=2,n ! Прямая прогонка
    T(i,2)=T(i,2)-T(i,1)*T(i-1,3)/T(i-1,2)
    T(i,4)=T(i,4)-T(i,1)*T(i-1,4)/T(i-1,2)
    T(i,1)=0d0
  end do
  T(n,4)=T(n,4)/T(n,2)
  do i=n-1,1,-1 ! Обратная прогонка
    T(i,4)=(T(i,4)-T(i,3)*T(i+1,4))/T(i,2)
    T(i,2)=1d0
    T(i,3)=0d0
  end do
  T(n,2)=1d0
!--------------------------------------------------------------------
return
end subroutine TDM

subroutine TDM2(n,T)
!--------------------------------------------------------------------
! Блок описания переменных
integer*4, intent(in)                 :: n ! Число строк
real*8, dimension(n,4), intent(inout) :: T ! Трёхдиагональная матрица
integer*4                             :: i ! Индекс текущей строки
!--------------------------------------------------------------------
! Основные циклы
  do i=2,n ! Прямая прогонка
    T(i,2)=T(i,2)-T(i,1)*T(i-1,3)/T(i-1,2)
    T(i,4)=T(i,4)-T(i,1)*T(i-1,4)/T(i-1,2)
    T(i,1)=0d0
  end do
  print *, T(:,4)
  pause
  T(n,4)=T(n,4)/T(n,2)
  do i=n-1,1,-1 ! Обратная прогонка
    T(i,4)=(T(i,4)-T(i,3)*T(i+1,4))/T(i,2)
    T(i,2)=1d0
    T(i,3)=0d0
  end do
  T(n,2)=1d0
!--------------------------------------------------------------------
return
end subroutine TDM2

subroutine k_omega_gamma(dx,k,n,i,y,V,U,rho,gam,gamit,kt,ktit,omega,omegait,mu,mut)
integer*4, intent(in)                :: k,i,n
real*8, dimension(k), intent(in)     :: mu
real*8, dimension(n,k), intent(in)   :: y, rho, V, U
real*8, dimension(n,k), intent(inout):: gam,gamit,kt,ktit,omega,omegait
real*8, dimension(2,k), intent(inout):: mut
integer*4                            :: j
real*8                               :: gammamax, ohm, rnu, rt, tomega, rc, fturb, egam, pgam, ggam, fgam
real*8                               :: sigmal,sigmagamma,dym,dyp,dx,dyp1,S,pk,sigmak,c_mu,sigmaomega
real*8                               :: comega1,comega2,dymk,Rel,Relit
real*8, dimension(k,4)               :: abcd
!--------------------------------------------------------------------
  c_mu=0.09d0
  comega1=5d0/9d0
  comega2=3d0/40d0
  sigmak=2d0
  sigmaomega=2d0
  sigmal=5d0
  sigmagamma=0.2d0
  Tu=0.065d0
  TI=100d0
  !Finite Difference Solution gam
do j=2,k-1
    dym=y(i,j)-y(i,j-1)
    dyp=y(i,j+1)-y(i,j)
    gammamax=1.1d0
    gam(i,j)=min(gam(i,j),1d0)
    mut(i,j)=rho(i,j)*kt(i,j)/omega(i,j)
    ohm=dabs(((dyp/dym*(U(i,j)-U(i,j-1))+dym/dyp*(U(i,j+1)-U(i,j)))/(dym+dyp)))
    rnu=rho(i,j)*(y(i,j))**(2d0)*ohm/(2.188d0*mu(j))
    rt=mut(i,j)/mu(j)
    tomega=rt*ohm/omega(i,j)
    rc=400d0-360d0*min(tomega/2d0,1d0)
    if ((rnu<=rc).or.(rnu>=100d0/0.7d0)) fgam=0d0
    if ((rnu>rc+4d0).and.(rnu<=100d0/0.7d0-1d0)) fgam=8d0
    fturb=exp(-(rnu*rt)**(1.2d0))
    if ((rnu<=18d0).or.(rnu>=100d0)) ggam=0d0
    if ((rnu>19d0).and.(rnu<=99d0)) ggam=7.5d0
    egam=ggam*fturb*ohm*(gam(i,j))**(1.5d0)
    pgam=fgam*ohm*(gammamax-gam(i,j))*sqrt(gam(i,j))
    abcd(j,1)=-dyp*rho(i,j)*V(i,j)/dym/(dyp+dym)-(mu(j-1)/sigmal+mut(i,j-1)/sigmagamma)/dym*2d0/(dyp+dym)
    abcd(j,2)=rho(i,j)*U(i,j)/dx+rho(i,j)*V(i,j)*(dyp-dym)/dyp/dym+2d0*(mu(j+1)/sigmal &
    +mut(i,j+1)/sigmagamma)/dyp/(dyp+dym)+2d0*(mu(j-1)/sigmal+mut(i,j-1)/sigmagamma)/dym/(dyp+dym)
    abcd(j,3)=rho(i,j)*V(i,j)*dym/dyp/(dyp+dym)-2d0*(mu(j+1)/sigmal+mut(i,j+1)/sigmagamma)/dyp/(dym+dyp)
    abcd(j,4)=rho(i,j)*U(i,j)*gam(i-1,j)/dx+rho(i,j)*(pgam-egam)

end do



dyp1=y(i,2)-y(i,1)
abcd(1,1)=0d0
abcd(1,2)=-1d0/dyp1
abcd(1,3)=1d0/dyp1
abcd(1,4)=0d0

abcd(k,1)=0d0
abcd(k,2)=1d0
abcd(k,3)=0d0
abcd(k,4)=1d0

call TDM(k,abcd)
Rel=0.5d0
Relit=1d0-Rel
gam(i,:)=Rel*abcd(:,4)+Relit*gamit(i,:)



!Finite Difference Solution kt
do j=2,k-1
    dym=y(i,j)-y(i,j-1)
    dyp=y(i,j+1)-y(i,j)
    mut(i,j)=rho(i,j)*kt(i,j)/omega(i,j)
    ohm=dabs(((dyp/dym*(U(i,j)-U(i,j-1))+dym/dyp*(U(i,j+1)-U(i,j)))/(dym+dyp)))
    S=sqrt(0.5)*((dyp/dym*(U(i,j)-U(i,j-1))+dym/dyp*(U(i,j+1)-U(i,j)))/(dym+dyp))
    pk=gam(i,j)*min(2d0*mut(i,j)*S**2d0/rho(i,j),kt(i,j))*abs(S/sqrt(3d0))
    abcd(j,1)=-dyp*rho(i,j)*V(i,j)/dym/(dyp+dym)-(mu(j-1)+mut(i,j-1)/sigmak)/dym*2d0/(dyp+dym)
    abcd(j,2)=rho(i,j)*U(i,j)/dx+rho(i,j)*V(i,j)*(dyp-dym)/dyp/dym+2d0*(mu(j+1)+mut(i,j+1)/sigmak)/dyp/(dyp+dym)&
    +2d0*(mu(j-1)+mut(i,j-1)/sigmak)/dym/(dyp+dym)+c_mu*rho(i,j)*omega(i,j)
    abcd(j,3)=rho(i,j)*V(i,j)*dym/dyp/(dyp+dym)-2d0*(mu(j)+mut(i,j)/sigmak)/dyp/(dym+dyp)
    abcd(j,4)=rho(i,j)*U(i,j)*kt(i-1,j)/dx+rho(i,j)*pk
end do

abcd(1,1)=0d0
abcd(1,2)=1d0
abcd(1,3)=0d0
abcd(1,4)=0d0

dymk=y(i,k)-y(i,k-1)
abcd(k,1)=-1d0/dymk
abcd(k,2)=1d0/dymk
abcd(k,3)=0d0
abcd(k,4)=0d0

call TDM(k,abcd)
Rel=0.05d0
Relit=1d0-Rel
kt(i,:)=Rel*abcd(:,4)+Relit*ktit(i,:)

!Finite Difference Solution omega
do j=2,k-1
    dym=y(i,j)-y(i,j-1)
    dyp=y(i,j+1)-y(i,j)
    mut(i,j)=rho(i,j)*kt(i,j)/omega(i,j)
    ohm=dabs(((dyp/dym*(U(i,j)-U(i,j-1))+dym/dyp*(U(i,j+1)-U(i,j)))/(dym+dyp)))
    S=sqrt(0.5d0)*((dyp/dym*(U(i,j)-U(i,j-1))+dym/dyp*(U(i,j+1)-U(i,j)))/(dym+dyp))
    abcd(j,1)=-dyp*rho(i,j)*V(i,j)/dym/(dyp+dym)-2d0*(mu(j-1)+mut(i,j-1)/sigmaomega)/dym/(dyp+dym)
    abcd(j,2)=rho(i,j)*U(i,j)/dx+rho(i,j)*V(i,j)*(dyp-dym)/dyp/dym+2d0*(mu(j+1)+mut(i,j+1)/sigmaomega)/dyp/(dyp+dym)&
    +2d0*(mu(j-1)+mut(i,j-1)/sigmaomega)/dym/(dyp+dym)+comega2*rho(i,j)*omega(i,j)
    abcd(j,3)=rho(i,j)*V(i,j)*dym/dyp/(dyp+dym)-2d0*(mu(j+1)+mut(i,j+1)/sigmaomega)/dyp/(dym+dyp)
    abcd(j,4)=rho(i,j)*U(i,j)*omega(i-1,j)/dx+2d0*rho(i,j)*comega1*S**2d0
end do

abcd(1,1)=0d0
abcd(1,2)=1d0
abcd(1,3)=0d0
abcd(1,4)=6d0*mu(2)/(comega2*rho(i,2)*(y(i,2))**(2d0))

dymk=y(i,k)-y(i,k-1)
abcd(k,1)=-1d0/dymk
abcd(k,2)=1d0/dymk
abcd(k,3)=0d0
abcd(k,4)=0d0

call TDM(k,abcd)
Rel=0.05d0
Relit=1d0-Rel
omega(i,:)=Rel*abcd(:,4)+Relit*omegait(i,:)

!print *, it,omega(i,2), gam(i,2), kt(i,2)
!if (it >=3999)    pause

mut(i,:)=rho(i,:)*kt(i,:)/omega(i,:)
return
end subroutine
