    subroutine Main(ip,aa,aw,s,dL,dC)
! EMSLN - Embedded Microstrip Coupled Lines, N-strip 
! ЗМПЛ - Заглубленные в диэлектрик связанные микрополосковые линии
    implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
    dimension betam(82),qwork(2004)
    dimension aa(6), aw(9), s(8)
    dimension z1(74), z2(74), z3(74), x3(74), x4(74),  x5(36), aC(289)
    dimension z1p(76),z2p(76),z3p(76),x3p(76),x4p(76),x5p(38),aCp(289)
    dimension aws(20),xs(20), bC(81), bCp(81), dC(81), dL(81)
    zinf = (1.d20,1.d20)
    zero = (0.d0, 0.d0)
    zi = (0.d0, 1.d0)
    pi = dacos(-1.d0)
    n  = aa(1)
    h1 = aa(2)
    h2 = aa(3)
    t  = aa(4)
    e1 = aa(5)
    e2 = aa(6)
    nn = 8*n+2
    do 10 i=1,n
        i8=8*(i-1);    sw=aw(i)/2
        do j=1,i-1; 
            sw=sw+aw(j)+s(j);
        enddo
        z1(i8+1) = sw;					betam(i8+1) = -0.5	! beta = alpha - 1
        z1(i8+2) = z1(i8+1) + zi*h1;	betam(i8+2) = -0.5
        z1(i8+3) = z1(i8+2) - aw(i)/2;  betam(i8+3) =  0.5
        z1(i8+4) = z1(i8+3) + zi*t;     betam(i8+4) =  0.5
        z1(i8+5) = z1(i8+4) + aw(i);	betam(i8+5) =  0.5
        z1(i8+6) = z1(i8+5) - zi*t;		betam(i8+6) =  0.5
        z1(i8+7) = z1(i8+2);			betam(i8+7) = -0.5
        z1(i8+8) = z1(i8+1);			betam(i8+8) = -0.5
10  continue
    sww=sw+aw(n)/2.
    sww2=sww/2.
    z1(8*n+1) = zinf
    z1(8*n+2) = -10000.
    if(n<4) z1(8*n+2) =-100.
    betam(8*n+1) = -2.
    betam(8*n+2) =  0.
    z10 = dcmplx(sww2, sww/8+(h1+t))
    nq = 13
    call qinit(nn,betam,nq,qwork)
    do k = 1,nn
        z2(k) = exp(dcmplx(0.d0, k-nn));
    enddo
    tol = 1.d-11
    iguess=1
    call scsolv(ip, iguess, tol, errest, nn, c, z2, z10, z1, betam,nq,qwork)
    z30 = zi;
    z301=-zi;
    do k = 1, nn-2
        z3(k)=(z2(k)*z301 - z30)/(z2(k) - 1.)
        x3(k)=dreal(z3(k));
    enddo
    do k = 1,nn-2
        x4(k) = 2 * (x3(k) - x3(1)) / (x3(nn-2) - x3(1)) - 1.;
    enddo
    do i=1,n
        i4=4*(i-1);			i8=8*(i-1);
        x5(i4+1)=x4(i8+1);	x5(i4+3)=x4(i8+7);
        x5(i4+2)=x4(i8+2);	x5(i4+4)=x4(i8+8);
    enddo
    M=1100
    n2=2*n-1
    call GHIONE(x5,aC,n2,M)
    call refor(aC,bC,n)
    call induc(bC,dL,n)
    nnp = 8*n+3
    aws1=aw(1)+s(1); aws2=aws1+aw(2); aws01=aw(1)/2.; aws02=aws1+aw(2)/2.
    zh1=zi*h1;       zh1t=zi*(h1+t);  zh1th2=zi*(h1+t+h2)
    sw = sww - aw(n)/2.
    do 20 i=1,n
        k=n+1-i;	i8=8*(i-1)
        z1p(i8+1) = sw + zh1t;				betam(i8+1) = .5
        z1p(i8+2) = z1p(i8+1) + aw(k)/2.;	betam(i8+2) = 1.5
        z1p(i8+3) = z1p(i8+2) - zi*t;		betam(i8+3) = 1.5
        z1p(i8+4) = z1p(i8+3) - aw(k);		betam(i8+4) = 1.5
        z1p(i8+5) = z1p(i8+4) + zi*t;		betam(i8+5) = 1.5
        z1p(i8+6) = z1p(i8+1);				betam(i8+6) = .5
        z1p(i8+7) = z1p(i8+1) + zi*h2;		betam(i8+7) = .5
        if(k/=1)then; 
            sw=sw-(aw(k)+aw(k-1))/2 -s(k-1)
            z1p(i8+8) = sw + zh1th2;			
            betam(i8+8) = .5;
        endif
20  continue
    z1p(nnp-3) =-zinf;					betam(nnp-3) = 0.
    z1p(nnp-2) = sww2;					betam(nnp-2) = 1.
    z1p(nnp-1) = zinf;					betam(nnp-1) = 0.
    z1p(nnp)   = z1p(7);				betam(nnp)   = .5
    betam=betam-1.
    z10p = dcmplx(sww2, 0.5 * h1  )
    nq = 13
    call qinit(nnp,betam,nq,qwork)
    do k = 1,nnp
        z2p(k) = exp(dcmplx(0.d0, k-nnp));
    enddo
    tol = 1.d-11
    iguess=1
    call scsolv(ip, iguess, tol, errest, nnp, c, z2p, z10p, z1p, betam,nq,qwork)
    z30p = zi;
    z301p=-zi;
    do k = 1, nnp-1
        z3p(k)=(z2p(k)*z301p - z30p)/(z2p(k) - 1.)
        x3p(k)=dreal(z3p(k));
    enddo
    do k = 1,nnp-1
        x4p(k) = 2 * (x3p(k) - x3p(1)) / (x3p(nnp-1) - x3p(1)) - 1.;
    enddo
    do i=1,2*n,2;
        j1 =(i-1)*4 + 1;
        j2 = j1 + 5;

        x5p(i) = x4p(j1);
        x5p(i+1) = x4p(j2);
    enddo
    x5p(2*n+1)=x4p(j2+2);	x5p(2*n+2)=x4p(j2+4)
    M=2000
    call GHIONE2(x5p,aCp,n,M)
    do i=1,n*n;
        bCp(i) = aCp(n*n+1-i);
    enddo
    dC = 8.854*( e2*bC + (e1-e2)*bCp)
    end