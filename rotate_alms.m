function alm_out=rotate_alms(lmax, alm, psi, theta, phi)
%function alm_out=rotate_alms(lmax, alm, psi,theta, phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converted from fortran to matlab 
%(from the healpix function "rotate_alm_d.f90") 
%J.Tolan (20141022)
%
%INPUTS:
%
%        lmax: maximum lmax to include in rotation. output will
%        contain entries only up to this lmax
%        alm: alm entries
%        psi, theta, phi: Euler angles
% 
%
%=========================================================
%Input: Complex array alm(p,l,m) with (l,m) in [0,lmax]^2, and p in [1,nd]
%Euler rotation angles psi, theta, phi in radians
%Output: Rotated array alm(p, l,m)
%
% Euler angle convention  is right handed, active rotation
% psi is the first rotation about the z-axis (vertical), in [-2pi,2pi]
% then theta about the ORIGINAL (unrotated) y-axis, in [-2pi,2pi]
% then phi  about the ORIGINAL (unrotated) z-axis (vertical), in [-2pi,2pi]
%
% Equivalently
% phi is the first rotation about the z-axis (vertical)
% then theta  about the NEW   y-axis (line of nodes)
% then psi    about the FINAL z-axis (figure axis)
% ---
% the recursion on the Wigner d matrix is inspired from the very stable
% double sided one described in Risbo (1996, J. of Geodesy, 70, 383)
% based on equation (4.4.1) in Edmonds (1957).
% the Risbo's scatter scheme has been repladed by a gather scheme for
% better computing efficiency
% the size of the matrix is divided by 2 using Edmonds Eq.(4.2.5) 
% to speed up calculations
% the loop on j has been unrolled for further speed-up
% EH, March--April 2005
%=========================================================


if (abs(psi) > 2.d0*pi || abs(phi) > 2.d0*pi || abs(theta) > 2.d0*pi)
  error('angles should be in radians')
end

nd = size(alm,1); 

na1=lmax+1;
na2=lmax+1;

if (na1 < (lmax+1) || na2 < (lmax+1))
  error('unconsistent alm array size and lmax')
end

d=zeros(2*lmax+2, lmax+2);
dd=zeros(2*lmax+2, lmax+2);
sqt=zeros(2*lmax+1,1);rsqrt=zeros(2*lmax+1,1);
alm1=zeros(nd, lmax+1);alm2=zeros(nd,lmax+1);


for i=1:2:lmax+1
  tsign(i)   =  1.0;
  tsign(i+1) = -1.0;
end
    
%initialization of square-root  table
    for i=1:2*lmax+1
       sqt(i) = sqrt(double(i-1));
    end

% initialisation of exponential table
for i=1:lmax+1
  exppsi(i)= complex(cos(psi*(i-1)), -sin(psi*(i-1)));
  expphi(i)= complex(cos(phi*(i-1)), -sin(phi*(i-1)));
end

% Note: theta has the correct sign.
p = sin(theta/2.d0);
q = cos(theta/2.d0);

disp('rotating_alms...')
for l=1:lmax+1
  ll=l+1; %index is incremented b one for -1 index in fortran
  ell=l-1; %ell value is one less than loop value
  
  % ------ build d-matrix of order l ------
  if (ell == 0)
    d(2,2) = 1.d0;
  end
  
  if(ell>=1)
  
  if (ell == 1)
    %initialize d-matrix degree 1/2
    dd(2,2)  =  q;
    dd(3,2)  = -p;
    dd(2,3)  =  p;
    dd(3,3)  =  q;
  end
  
  if (ell >=2)
  %  l - 1 --> l - 1/2
  jj = 2*ell-1;
  j=jj+1;
  rsqt = flipud(sqt(1:j));
  fj = jj;
  qj = q / fj;
  pj = p / fj;

  for kk = 0:jj/2 % keep only m' <= -1/2
    k=kk+1;
    dd(2:j+1,k+1) = rsqt(1:j).* ( d(2:j+1,k+1).* (sqt(jj-kk+1)  * qj)...
                          +       d(2:j+1,k).*   (sqt(kk+1) *    pj) )...
                 +  sqt(1:j).*  ( d(1:j,k).*     (sqt(kk+1) *    qj)... 
	                      -   d(1:j,k+1).*   (sqt(jj-kk+1)  * pj) );
  end % loop on k
  
  
  % l=half-integer, reconstruct m'= 1/2 by symmetry
  hhj = ell-1;
  hj=hhj+1;
  if (mod(ell,2) == 0) 
    for kk = 0:2:jj-1
      k=kk+1;
      dd(k+1,   ll) =   dd(j-k+2,   hj+1);
      dd(k+2, ll) = - dd(j-k+1, hj+1);
    end
  else
    for k = 1:2:j-1
      dd(k+1,   ll) = - dd(j-k+2,   hj+1);
      dd(k+2, ll) =   dd(j-k+1, hj+1);
    end
  end
  
  end %if ell>=2

  %  l - 1/2 --> l
  jj = 2*ell;
  j=jj+1;
  rsqt = flipud(sqt(1:j));
  fj = jj;
  qj = q / fj;
  pj = p / fj;
  
  for kk = 0:j/2 % keep only m' <= 0
    k=kk+1;
    d (2:j+1,k+1) = rsqt(1:j).*  ( dd(2:j+1,k+1).*  (sqt(jj-kk+1)  * qj)...
                               +   dd(2:j+1,k).*    (sqt(kk+1)   *  pj) )...
               + sqt(1:j).*      ( dd(1:j,k).*      (sqt(kk+1)   *   qj)... 
	                         - dd(1:j,k+1).*    (sqt(jj-kk+1)  *  pj) );
  end % loop on k

  end %if ell>=1
  

  % ------- apply rotation matrix -------
 % kd=1;
 for kd=1:nd
   alm_s=squeeze(alm(kd,:,:));
   alm1(kd,1:l)  = alm_s(l,1:l).*exppsi(1:l);
 end
 
  % m = 0
  for kd=1:nd
    alm2(kd,1:l) = (alm1(kd,1).*d(ell+2:2*ell+2,l+1)).';
  end
  
  if(ell~=0)
    for mm = 1:l
      for kd=1:nd
	alm2(kd, mm) = alm2(kd,mm) + sum(alm1(kd,2:l).* flipud(d(2:ell+1,l-mm+2)).')...
	    + conj(sum(alm1(kd,2:l).* (tsign(2:l).* d(ell+3:2*ell+2,l-mm+2).')));
      end
    end
  end  

  %add to the alm entries 
  for kd=1:nd
    alm_out(kd,l,1:l)=alm2(kd,1:l).*expphi(1:l);
  end
  
    %clear these in this ell loop  
    clear alm1;
    clear alm2;
       

end % loop on ll
  
disp('done rotating_alms')
return

