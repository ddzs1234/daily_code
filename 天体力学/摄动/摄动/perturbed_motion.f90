!author Yang Shihai & Zhai Sai
!本程序先通过线性长期摄动方程解出轨道根数的变化，然后再用一般的n体问题解来验证线性长期摄动方程解的准确性
!程序输出为线性撑起摄动方程解的轨道根数变化以及一般n体问题解的位置速度及其轨道根数变化，然后用python程序进行绘图
    program main
    
    implicit none
    
    !variables
    integer(kind=4),parameter::num=2 !考虑行星的个数
    integer(kind=4) i,j,k,count1,count2,num1,n,targ,count0  !targ为选取的要保存数据的行星，targ=1时为木星，targ=2时为土星
    real(kind=8) b(13,13),b0(13,13),b1(13,13),a(13,1),a0(13,1),a1(13,1),c(13,1),c0(13,1),c1(13,1),c2(13,1),c10(13,1),c11(13,1)
    real(kind=8) t,x(6*num,1),vector(4*num,1),tend,h,err,err1,vector0(4*num,1)   !err是线性长期摄动积分过程的允许误差，err1是用普通的牛顿运动方程的允许误差
    real(kind=8) ::pi=acos(-1d0),g_const=6.67428d-11,au=149597870691d0,m0=1.989d30
    real(kind=8) eccen(num,1),inclination(num,1),omega(num,1),pomega(num,1),lamda(num,1),semi_axis(num,1),mass(num,1),n0(num,1)
    real(kind=8) alpha(num,num),a_mat(num,num),b_mat(num,num),tran_mat(4*num,4*num),c00(num,num),c_i(num,num)
    real(kind=8) a1_mat(num,num),i3(3,3),e0(num,1),mtemp(num,1),r,position(3,2),velocity(3,2),elements1(2,2)
    real(kind=8) r0,v0_2,semi,n2,hx,hy,hz,h_0,cos_e,sin_e,temp1,temp0,sign1,x0(6*num,1)
    real(kind=8) vector1(4,1),elements(4,2),temp(3,1)
    real(kind=8),external::laplace_coef
    common /group1/b,a,c,c1
    common /group2/tran_mat
    common /group3/i3,g_const,mass,m0,num1

    !给一些变量赋值
    data b0 /0.0,2.0,1.0,1.0,5.0,1.0,-25.0,31.0,2.0,-91.0,2383.0,3.0,-1777.0,&
        0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,-25.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
        0.0,0.0,0.0,0.0,0.0,25.0,1.0,125.0,0.0,-53.0,23.0,-341.0,0.0,-341.0,0.0,0.0,0.0,0.0,0.0,&
        1.0,-65.0,61.0,704.0,-976.0,4496.0,0.0,4496.0,0.0,0.0,0.0,0.0,0.0,0.0,125.0,-2.0,&
        -107.0,311.0,-301.0,-6.0,-289.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,13.0,67.0,-19.0,2133.0,&
        -3.0,2193.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,17.0,45.0,-3.0,51.0,0.0,0.0,0.0,0.0,0.0,0.0,&
        0.0,0.0,0.0,-1.0,45.0,3.0,33.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,18.0,6.0,12.0,&
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,&
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
    data b1 /1.0,27.0,36.0,24.0,12.0,20.0,108.0,300.0,1.0,108.0,4100.0,205.0,4100.0,&
        1.0,1.0,12.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,8.0,16.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,&
        1.0,1.0,1.0,1.0,1.0,16.0,4.0,108.0,1.0,6.0,108.0,164.0,1.0,164.0,1.0,1.0,1.0,1.0,1.0,&
        5.0,27.0,225.0,45.0,135.0,1025.0,1.0,1025.0,1.0,1.0,1.0,1.0,1.0,1.0,54.0,9.0,&
        9.0,54.0,82.0,41.0,82.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,900.0,90.0,60.0,4100.0,&
        205.0,4100.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,6.0,82.0,41.0,82.0,1.0,1.0,1.0,1.0,1.0,1.0,&
        1.0,1.0,1.0,12.0,162.0,41.0,164.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,41.0,41.0,41.0,&
        1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,&
        1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
	data a0 /0.0,2.0,1.0,1.0,5.0,1.0,5.0,1.0,2.0,1.0,1.0,0.0,1.0/
    data a1 /1.0,27.0,9.0,6.0,12.0,2.0,6.0,6.0,3.0,3.0,1.0,1.0,1.0/
    data c0 /41.0,0.0,0.0,0.0,0.0,34.0,9.0,9.0,9.0,9.0,41.0,0.0,0.0/
    data c2 /840.0,1.0,1.0,1.0,1.0,105.0,35.0,35.0,280.0,280.0,840.0,1.0,1.0/
    data c10 /0.0,0.0,0.0,0.0,0.0,34.0,9.0,9.0,9.0,9.0,0.0,41.0,41.0/
    data c11 /1.0,1.0,1.0,1.0,1.0,105.0,35.0,35.0,280.0,280.0,1.0,840.0,840.0/
    data eccen /0.04839266,0.05415060/
    data inclination /1.30530,2.48446/
    data omega /100.55615,113.71504/
    data pomega /14.75385,92.43194/
    data lamda /34.40438,49.94432/
    data semi_axis /5.20336301,9.53707032/
    data mass /18986d23,5684.6d23/
    data i3 /1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/
    
    num1=num
    tend=100000.0*365.0*86400.0  !结束时间
    b=b0/b1
    b=transpose(b)
    a=a0/a1
    c=c0/c2
    c1=c10/c11
    inclination=inclination*pi/180.0
    omega=omega*pi/180.0
    pomega=pomega*pi/180.0
    lamda=lamda*pi/180.0
    semi_axis=semi_axis*au
    n0=(g_const*(m0+mass)/semi_axis**3)**0.5
    !得到长期线性摄动方程中的状态转移矩阵
    alpha=0d0
    do i=1,num
        do j=i+1,num
            alpha(i,j)=semi_axis(i,1)/semi_axis(j,1)
        end do
    end do
    alpha=alpha+transpose(alpha)
    c00=0d0
    do i=1,num
        do j=i+1,num
            c00(i,j)=laplace_coef(2d0,1.5d0,alpha(i,j))/laplace_coef(1.0d0,1.5d0,alpha(i,j))
        end do
    end do
    c00=c00+transpose(c00)
    c_i=0d0
    do i=1,num
        do j=1,num
            if (i<j) then
                c_i(i,j)=0.25*n0(i,1)*mass(j,1)/(m0+mass(i,1))*alpha(i,j)**2*laplace_coef(1d0,1.5d0,alpha(i,j))
            else if (i>j) then
                c_i(i,j)=0.25*n0(i,1)*mass(j,1)/(m0+mass(i,1))*alpha(i,j)*laplace_coef(1d0,1.5d0,alpha(i,j))
            end if
        end do
    end do
    a1_mat=0d0
    do i=1,num
        a1_mat(i,i)=sum(c_i(i,:))
    end do
    a_mat=a1_mat-c00*c_i
    b_mat=c_i-a1_mat
    tran_mat=0d0
    tran_mat(1:num,num+1:2*num)=a_mat
    tran_mat(num+1:2*num,1:num)=-a_mat
    tran_mat(2*num+1:3*num,3*num+1:4*num)=b_mat
    tran_mat(3*num+1:4*num,2*num+1:3*num)=-b_mat
    
    !设定长期线性摄动方程所用初始值
    vector=0d0
    vector(1:num,1:1)=eccen*sin(pomega)
    vector(num+1:2*num,1:1)=eccen*cos(pomega)
    vector(2*num+1:3*num,1:1)=inclination*sin(omega)
    vector(3*num+1:4*num,1:1)=inclination*cos(omega)
    count1=1
    n=4*num
    err=1e-16
    h=8640000.0
    t=0.0
    vector0=vector
    open(unit=14,file="elements.txt")
    write(14,1012) t/365.0/86400.0,eccen(1,1),pomega(1,1)*180.0/pi,inclination(1,1)*180.0/pi,omega(1,1)*180.0/pi,&
        eccen(2,1),pomega(2,1)*180.0/pi,inclination(2,1)*180.0/pi,omega(2,1)*180.0/pi
    do while (t<tend) !计算t=0以后的轨道根数变化
        call RKF78(h,t,vector,err,n)
        t=t+h
        do targ=1,2
            vector1(1:4,1)=vector(targ::num,1)
            !接下来通过得到的考虑奇点的根数来计算e,i,pomega以及升交点经度
            elements(1,targ)=(vector1(1,1)**2+vector1(2,1)**2)**0.5
            temp0=acos(vector1(2,1)/elements(1,targ))
            call sgn(sign1,vector1(1,1))
            elements(2,targ)=(temp0+(pi-temp0)*sign1*(sign1-1d0))*180.0/pi
            elements(3,targ)=(vector1(3,1)**2+vector1(4,1)**2)**0.5
            temp0=acos(vector1(4,1)/elements(3,targ))
            call sgn(sign1,vector1(3,1))
            elements(4,targ)=sign1*temp0*180.0/pi
            elements(3,targ)=elements(3,targ)*180.0/pi
        end do
        write(14,1012) t/365.0/86400.0,elements(1:4,1),elements(1:4,2)
    end do
    h=-h
    t=0.0
    vector=vector0
    do while (t>-tend) !计算t=0以前的轨道根数变化
        call RKF78(h,t,vector,err,n)
        t=t+h
        do targ=1,2
            vector1(1:4,1)=vector(targ::num,1)
            !接下来通过得到的考虑奇点的根数来计算e,i,pomega以及升交点经度
            elements(1,targ)=(vector1(1,1)**2+vector1(2,1)**2)**0.5
            temp0=acos(vector1(2,1)/elements(1,targ))
            call sgn(sign1,vector1(1,1))
            elements(2,targ)=(temp0+(pi-temp0)*sign1*(sign1-1d0))*180.0/pi
            elements(3,targ)=(vector1(3,1)**2+vector1(4,1)**2)**0.5
            temp0=acos(vector1(4,1)/elements(3,targ))
            call sgn(sign1,vector1(3,1))
            elements(4,targ)=sign1*temp0*180.0/pi
            elements(3,targ)=elements(3,targ)*180.0/pi
        end do
        write(14,1012) t/365.0/86400.0,elements(1:4,1),elements(1:4,2)
    end do
    close(14)
    print *, "摄动方程求解结束"
    !摄动方程求解过程结束
    !用牛顿运动方程求解过程开始
    mtemp=lamda-pomega
    do i=1,num
        if (mtemp(i,1)<0) then
            mtemp(i,1)=mtemp(i,1)+2d0*pi
        end if
    end do
    call m2e(e0,eccen,mtemp,num)
    do i=1,num        !通过轨道根数的资料转换成初始坐标
        temp(1,1)=semi_axis(i,1)*(cos(e0(i,1))-eccen(i,1))
        temp(2,1)=semi_axis(i,1)*(1d0-eccen(i,1)**2)**0.5*sin(e0(i,1))
        temp(3,1)=0d0
        call rz(omega(i,1)-pomega(i,1),temp)
        call rx(-inclination(i,1),temp)
        call rz(-omega(i,1),temp)
        x(3*i-2:3*i,1:1)=temp
        r=dot_product(temp(:,1),temp(:,1))**0.5
        temp(1,1)=-semi_axis(i,1)**2*n0(i,1)/r*sin(e0(i,1))
        temp(2,1)=semi_axis(i,1)**2*n0(i,1)/r*(1d0-eccen(i,1)**2)**0.5*cos(e0(i,1))
        temp(3,1)=0d0
        call rz(omega(i,1)-pomega(i,1),temp)
        call rx(-inclination(i,1),temp)
        call rz(-omega(i,1),temp)
        x(3*num+3*i-2:3*num+3*i,1:1)=temp
    end do
    t=0.0
    n=6*num
    h=8640000.0
    err1=1d-14*au
    x0=x
    open(unit=10,file="position.txt")
    open(unit=12,file="velocity.txt")
    open(unit=16,file="elements1.txt")
    do targ=1,2
        position(:,targ:targ)=x(3*targ-2:3*targ,1:1)
        velocity(:,targ:targ)=x(3*num+3*targ-2:3*num+3*targ,1:1)
    end do
    write(10,1010) t/365.0/86400.0,position
    write(12,1010) t/365.0/86400.0,velocity
    write(16,1012) t/365.0/86400.0,eccen(1,1),inclination(1,1)*180.0/pi,eccen(2,1),inclination(2,1)*180.0/pi
    do while (t<tend)  !计算t=0以后的轨道根数变化
        call RKF78_1(h,t,x,err1,n)
        t=t+h
        print *, t/365.0/86400.0
        do targ=1,2
            position(1:3,targ)=x(3*targ-2:3*targ,1)
            velocity(1:3,targ)=x(3*num+3*targ-2:3*num+3*targ,1)
        end do
        write(10,1010) t/365.0/86400.0,position
        write(12,1010) t/365.0/86400.0,velocity
        !通过计算的坐标数据得到轨道根数
        do targ=1,2
            r0=(position(1,targ)**2+position(2,targ)**2+position(3,targ)**2)**0.5
            v0_2=velocity(1,targ)**2+velocity(2,targ)**2+velocity(3,targ)**2
            semi=1.0/(2.0/r0-v0_2/(g_const*(m0+mass(targ,1))))
            h_0=((position(2,targ)*velocity(3,targ)-position(3,targ)*velocity(2,targ))**2+&
               (position(3,targ)*velocity(1,targ)-position(1,targ)*velocity(3,targ))**2+&
               (position(1,targ)*velocity(2,targ)-position(2,targ)*velocity(1,targ))**2)**0.5
            hz=position(1,targ)*velocity(2,targ)-position(2,targ)*velocity(1,targ)
            elements1(1,targ)=(1d0-h_0**2/(g_const*(m0+mass(targ,1))*semi))**0.5
            elements1(2,targ)=acos(hz/h_0)*180.0/pi
        end do
        write(16,1012) t/365.0/86400.0,elements1
    end do
    t=0.0
    h=-h
    x=x0
    do while (t>-tend)  !计算t=0以前的轨道根数变化
        call RKF78_1(h,t,x,err1,n)
        t=t+h
        print *, t/365.0/86400.0
        do targ=1,2
            position(1:3,targ)=x(3*targ-2:3*targ,1)
            velocity(1:3,targ)=x(3*num+3*targ-2:3*num+3*targ,1)
        end do
        write(10,1010) t/365.0/86400.0,position
        write(12,1010) t/365.0/86400.0,velocity
        !通过计算的坐标数据得到轨道根数
        do targ=1,2
            r0=(position(1,targ)**2+position(2,targ)**2+position(3,targ)**2)**0.5
            v0_2=velocity(1,targ)**2+velocity(2,targ)**2+velocity(3,targ)**2
            semi=1.0/(2.0/r0-v0_2/(g_const*(m0+mass(targ,1))))
            h_0=((position(2,targ)*velocity(3,targ)-position(3,targ)*velocity(2,targ))**2+&
               (position(3,targ)*velocity(1,targ)-position(1,targ)*velocity(3,targ))**2+&
               (position(1,targ)*velocity(2,targ)-position(2,targ)*velocity(1,targ))**2)**0.5
            hz=position(1,targ)*velocity(2,targ)-position(2,targ)*velocity(1,targ)
            elements1(1,targ)=(1d0-h_0**2/(g_const*(m0+mass(targ,1))*semi))**0.5
            elements1(2,targ)=acos(hz/h_0)*180.0/pi
        end do
        write(16,1012) t/365.0/86400.0,elements1
    end do   
    close(10)
    close(12)
    close(16)
1010 format(e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10)
1012 format(e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10,1x,e17.10)
    stop
    end program main
    
    function laplace_coef(j,s,alpha)  !计算拉普拉斯系数的函数
    
    implicit none
    real(kind=8)::laplace_coef,f1,func,j,s,pi=acos(-1d0),alpha
    real(kind=8),external::integration
    real(kind=8)::x,a,b
    laplace_coef=integration(0d0,2.0*pi,j,s,alpha,1d-10)/pi
    return
    end
    
    
    recursive function integration(x,y,j,s,alpha,err)    !自适应积分函数，通过递归得到最后积分的精确值
    implicit none
    real(kind=8)::f1,func,j,s,integration,x0,a,b,x,y,err,alpha
    f1(x0)=cos(j*x0)/(1.0-2.0*alpha*cos(x0)+alpha**2)**s
    func(a,b)=(b-a)/6.0*(f1(a)+4.0*f1((a+b)/2.0)+f1(b))
    if (abs(func(x,y)-func(x,(x+y)/2)-func((x+y)/2,y))<err) then
        integration=func(x,y)
    else
        integration=integration(x,(x+y)/2,j,s,alpha,err/2)+integration((x+y)/2,y,j,s,alpha,err/2)
    end if
    return
    end
    
    
    subroutine sgn(sign0,x)     !符号函数
    
    implicit none
    real(kind=8) sign0,x
    if (x<0) then
        sign0=-1d0
    else if (x>0) then
        sign0=1d0
    else
        sign0=0d0
    end if
    return
    end

    
    subroutine right_func(f,t,x0,n)    !长期线性摄动方程的右函数
    
    implicit none
    integer n
    integer,parameter:: num1=2
    real(kind=8)::f(n,1),x0(n,1),tran_mat(4*num1,4*num1),t
    common /group2/tran_mat
    
    f=matmul(tran_mat,x0)
    return
    end
    
    
    recursive subroutine RKF78(h,t,x,err,n)      !长期线性摄动方程的RKF78子程序
    
    implicit none
    real(kind=8) h,t,err,ttemp,h1,t1
    integer i,n
    real(kind=8) x(n,1),f(n,1),b(13,13),c(13,1),c1(13,1),a(13,1),k(n,13),x0(n,1),x1(n,1),temp(n,1),xtemp(n,1),temp_delta0
	common /group1/b,a,c,c1
    k=0
    do i=1,13
        xtemp=x+h*matmul(k,b(:,i:i))
        ttemp=t+a(i,1)*h
        call right_func(f,ttemp,xtemp,n)
        k(:,i:i)=f
	end do
	x0=x+h*matmul(k,c)
    x1=x+h*matmul(k,c1)
    temp=x1-x0
    temp_delta0=dot_product(temp(:,1),temp(:,1))**0.5/abs(h)
    if (temp_delta0>err) then
        h1=h/2
        call RKF78(h1,t,x,err,n)
        t1=t+h1
        call RKF78(h1,t1,x,err,n)
	else
        x=x0
    end if
	return
    end
    
    
    subroutine right_func_1(f,t,x0,n)     !一般的n体问题解法方程右函数
    
    implicit none
    integer n,num
    integer,parameter::num1=2
    real(kind=8)::f(n,1),x0(n,1),mat0(n,n),r_3(num,1),rij_3(num,num),mattemp(num,1),g_const,i3(3,3),t,mass(num1,1),m0
    integer i,j,k
    common /group3/i3,g_const,mass,m0,num
    
    mat0=0d0
    do i=1,num
        mat0(3*i-2:3*i,3*num+3*i-2:3*num+3*i)=i3
        r_3(i,1)=dot_product(x0(3*i-2:3*i,1),x0(3*i-2:3*i,1))**1.5
    end do
    rij_3=0d0
    do i=1,num
        do j=1,i-1
            rij_3(i,j)=dot_product(x0(3*i-2:3*i,1)-x0(3*j-2:3*j,1),x0(3*i-2:3*i,1)-x0(3*j-2:3*j,1))**1.5
        end do
    end do
    rij_3=rij_3+transpose(rij_3)
    do i=1,num
        do j=1,num
            if (i==j) then
                mattemp=0d0
                do k=1,num
                    if (k/=i) then
                        mattemp(k,1)=mass(k,1)/rij_3(i,k)
                    end if
                end do
                mat0(3*num+3*i-2:3*num+3*i,3*i-2:3*i)=(-g_const*(m0+mass(i,1))/r_3(i,1)-g_const*sum(mattemp))*i3
            else
                mat0(3*num+3*i-2:3*num+3*i,3*j-2:3*j)=g_const*mass(j,1)*(1.0/rij_3(i,j)-1.0/r_3(j,1))*i3
            end if
        end do
    end do
    
    f=matmul(mat0,x0)
    return
    end
    
    
    recursive subroutine RKF78_1(h,t,x,err,n)   !一般的n体问题方程的RKF78解法
    
    implicit none
    real(kind=8) h,t,err,ttemp,h1,t1
    integer i,n
    real(kind=8) x(n,1),f(n,1),b(13,13),c(13,1),c1(13,1),a(13,1),k(n,13),x0(n,1),x1(n,1),temp(n,1),xtemp(n,1),temp_delta0
	common /group1/b,a,c,c1
    k=0
    do i=1,13
        xtemp=x+h*matmul(k,b(:,i:i))
        ttemp=t+a(i,1)*h
        call right_func_1(f,ttemp,xtemp,n)
        k(:,i:i)=f
	end do
	x0=x+h*matmul(k,c)
    x1=x+h*matmul(k,c1)
    temp=x1-x0
    temp_delta0=dot_product(temp(:n/2,1),temp(:n/2,1))**0.5/abs(h)
    if (temp_delta0>err) then
        h1=h/2
        call RKF78_1(h1,t,x,err,n)
        t1=t+h1
        call RKF78_1(h1,t1,x,err,n)
	else
        x=x0
    end if
	return
    end
    
    
    subroutine m2e(e,eccen,m,n)    !将平近点角转化成偏近点角
    
    implicit none
    integer(kind=4) n,i
    real(kind=8)::e(n,1),eccen(n,1),m(n,1),sigma=1d-30,e0,e1,delta
    do i=1,n
        delta=1d0
        e0=0d0
        e1=0d0
        do while (delta>sigma)
            e1=e0-(e0-eccen(i,1)*sin(e0)-m(i,1))/(1-eccen(i,1)*cos(e0))
            delta=abs(e1-e0)
            e0=e1
        end do
        e(i,1)=e1
    end do
    return
    end
    
    
    subroutine rx(theta,x0)    !绕x轴旋转
    real(kind=8)::theta,x0(3,1),rot(3,3)
    rot=0d0
    rot(1,1)=1d0
    rot(2,2)=cos(theta)
    rot(3,3)=rot(2,2)
    rot(2,3)=sin(theta)
    rot(3,2)=-rot(2,3)
    x0=matmul(rot,x0)
    return
    end
    
    
    subroutine ry(theta,x0)     !绕y轴旋转
    real(kind=8)::theta,x0(3,1),rot(3,3)
    rot=0d0
    rot(2,2)=1d0
    rot(1,1)=cos(theta)
    rot(3,3)=rot(1,1)
    rot(3,1)=sin(theta)
    rot(1,3)=-rot(3,1)
    x0=matmul(rot,x0)
    return
    end
    
    
    subroutine rz(theta,x0)       !绕z轴旋转
    real(kind=8)::theta,x0(3,1),rot(3,3)
    rot=0d0
    rot(3,3)=1d0
    rot(1,1)=cos(theta)
    rot(2,2)=rot(1,1)
    rot(1,2)=sin(theta)
    rot(2,1)=-rot(1,2)
    x0=matmul(rot,x0)
    return
    end
 
    
