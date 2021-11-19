program main
  implicit none
  ! Блок описания переменных
real*8,   parameter      :: alpha=22.6d-6   ! коэффициент линейного теплового расширения 1/град
real*8,   parameter      :: den=16630d0
real*8,   parameter      :: hf=12.6d3 ! Дж/моль
real*8,   parameter      :: cp=250d0       ! Дж/кг*град
real*8                   :: dt,dr,dts       ! Шаг сетки по осям x,y
integer*4,parameter      :: n=9000,k=1400 ! Число узлов по осям x,y
!real*8,   dimension(n,1:k) :: T,r         ! Поля скоростей u,v
real*8,   dimension(n)   :: s
real*8,   dimension(k)   :: timeMelt, delCron
integer*4, dimension(k)  :: keyMelt
integer*4                :: i,j,js,bc,l,y,ms,g  ! Индексы текущего узла
real*8,   dimension(k,4) :: abcd        ! Трехдиагональная матрица
real*8,   dimension(7)   :: times
real*8,   dimension(7,k) :: Temps,eta
!--------------------------------------------------------------------
! Блок итерационных
character (60)::filenameProfiles
real*8,   dimension(n,k) :: Tit,T        ! Поля итерационных
real*8,   dimension(n,0:k) :: r                                            ! переменных U,V
real*8                   :: epsitT, epsidt         ! Погрешность определения
                                           ! скоростей U,V
real*8,   parameter      :: eps = 1d-10     ! Максимальная погрешность
integer*4                :: it             ! Номер итерации
real*8                   :: maxeps         ! Функция расчёта
                                           ! максимальной погрешности
character (60)::filenameProfiles2
!--------------------------------------------------------------------
! Блок расчета параметров сетки
  real*8 :: Rc, Tsol,Time,q,kt, Tm, Tout, delta, qliq, a,b,dtit,Timei,Timewrite,location
  real*8 :: Qz, St,zetast,tetast, Qov,flux

  open (55,file='timemf2.9.txt')
  write (55,91) 'tau'
  91 format (1(2x,A16))

  open (11,file='s(time)mf2.9.txt')
  write (11,91)'zeta'
  90 format(1(2x,A16))

  !open (56,file='temps1.txt')
  !write (56,91) 'temp'
  !91 format (1(2x,A16))

  !Important Parameters
  St=1d0 !cp*(Tm-Tout)/hf
  flux=2.9d0

  ms=1 !1 for melting, 2 for solidification
  bc=2 !1 for fixed temperature, 2 for fixed flux

  y=1 !Book keeping, sees what the first calculated dt is
  l=0 !sets js later
  g=0
  Tout=0d0 !1300d0
  Time=1d0 ! Seconds
  dt=1d-6  !Time/(n-1)
  dr=1d0/k!1d0/(k-1)!Rc/(k-1)
  if (bc==1)then
    Qz=7d0
  else
    Qz=9d0
  end if
  if(bc==1)then
    Tm=1d0
  else
    Tm=0d0
  end if
  if(bc==2)then
    St=1d0
  end if
!--------------------------------------------------------------------
! Блок начальных данных и граничных условий
  T=0d0         ! Обнуление массива для T
  r=0d0
  timeMelt(:)=0d0
  keyMelt(:)=0
  Timei=0d0         ! Обнуление массива для V
  !T(1,:)=0d0    ! 0d0 for melting 1d0 for solidification
  location=0d0
  do j=1,k !Set initial condition
    T(1,j)=-9d0/8d0*(location)**2d0!+9d0/6d0*(0.99d0)
    location=location+dr
  end do

  !do j=2,k-1
  !  T(1,j)=1d0/6d0*0.99d0**2d0*Qz-1d0/6d0*location**2d0*Qz+1d0
  !  location=location+dr
  !end do
  !T(:,k)=Tout
  !I need to make sure the initial state satisfies the BC
  !How should I account for this here when it's a flux condition?
  T(1,k)=-dr*flux+T(1,k-1)          ! T на внешней стенке цилиндра
!--------------------------------------------------------------------
! Блок начальных данных для итерационных переменных
  Tit = 0d0
  dtit= 0d0
  it = 0

  Timewrite=2.305d0 !0.2964d0+30d0 !0.3d0
!--------------------------------------------------------------------
! Основной цикл расчета
  if(ms==1)then
    js=0
  else
    js=k
  end if
  do i=2,n            ! Перебор по x
    T(i,:) = T(i-1,:) ! Первое приближение по U
    ! ---------------------------------------------------------------
    ! Итерационный цикл

     Do j=2,k-1
       r(i,j)=r(i,j-1)+dr
       if(ms==1)then
          if (js>0) js=js+1
       else
          if (js<k) js=js-1
       end if
       if ((ms==1).and.(js>0)) exit
       if ((ms==2).and.(js<k)) exit
       if (ms==1) then
        if ((T(i,j)<=Tm).and.(T(i,j-1)>Tm) ) then
          js=j
        end if
       else
        if ((T(i,j)>=Tm).and.(T(i,j+1)<Tm) ) then
          js=j
        end if
       end if
     end do

    do it =1,400
      Tit(i,:) = T(i,:) ! Обновление итерационных переменных по T
      ! ---------------------------------------------------------------
      ! Расчёт T

      if (ms==1) l=0
      if (ms==2) l=k
      If (js==(l)) then
         dt=dt !Time/(n-1)
      else
         a=1d0/(dr+dr**2d0/r(i,js))
         b=1d0/(dr-dr**2d0/r(i,js))
         dts=(0.5d0*dr**2d0*(a+b)*(Tm-T(i-1,js))+dr/St)/(b*(T(i,js+1)-Tm)-a*(Tm-T(i,js-1))+Qz*0.5d0*dr**2d0*(a+b))
         if(y==1)then
            dts=dabs(dts)
         end if
         y=y+1
         if (ms==1) dt=dts
         if (ms==2) dt=dabs(dts)
      end if

      do j=2,k-1       ! Перебор по y
        r(i,j)=r(i,j-1)+dr
        abcd(j,1)= 1d0/dr**2d0-1d0/r(i,j)/dr
        abcd(j,2)= -2d0/dr**2d0-1d0/dt
        abcd(j,3)= 1d0/dr**2d0+1d0/r(i,j)/dr
        abcd(j,4)=  -T(i-1,j)/dt -Qz
      end do

      r(i,k)=r(i,k-1)+dr

      ! В центре цилиндра j=1
      abcd(1,1)=0d0
      abcd(1,2)=-1d0/dr
      abcd(1,3)=1d0/dr
      abcd(1,4)=0d0

      abcd(js,1)=0d0
      abcd(js,2)=1d0
      abcd(js,3)=0d0
      abcd(js,4)=Tm
  if(bc==1)then
      abcd(k,1)=0d0
      abcd(k,2)=1d0
      abcd(k,3)=0d0
      abcd(k,4)=T(i,k)
  end if

  if(bc==2)then
      abcd(k,1)=-1d0/dr
      abcd(k,2)=1d0/dr
      abcd(k,3)=0d0
      abcd(k,4)=-flux
  end if


      call TDM(k,abcd) ! Вызов подпрограммы расчёта
                      ! трёхдиагональной матрицы
      T(i,:)=0.5d0*abcd(:,4)+0.5d0*Tit(i,:)
      ! ---------------------------------------------------------------
      epsidt=dabs(dt-dtit)
      epsitT = maxeps(T,Tit,n,k,i)
      dtit=dt
        if ((epsitT <= eps)) exit
    end do

    Timei=Timei+dt
    print *, i, it, 'js=', js, dt, Timei, r(i,js)
    ! ---------------------------------------------------------------
    if (ms==1)then
        if ((it>=400) .or. (dt<0d0)) Print*, 'it=401'
        if ((it>=400) .or. (dt<0d0)) exit
    else
        if ((it>=400) .or. (dts>0d0)) Print*, 'it=401'
        if ((it>=400) .or. (dts>0d0)) exit
    end if
    write (11,101) r(i,js)
    101 format(3(2x,E16.8))
    Write (55,95) Timei
    95 format (1(5x,E16.8))
    !if((Timei>=0.26233994E+00).and.(g==0))then
    !    g=1
    !    pause
    !    write(56,95) T(i,:)
    !    95 format (1(5x,E16.8))
    !end if

    zetast=r(i,js)
    tetast=T(i,1)
  end do
!--------------------------------------------------------------------

   Print*, 'Q=', Qz, 'St=', St, zetast,tetast
   100 format(3(35x,E16.8))
!--------------------------------------------------------------------
end program main

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
!  if (step==10) Print*, T(:,:)
!--------------------------------------------------------------------
return
end subroutine TDM
!  re_i = system("pause")
!end
