; routines for calculation of hydrogenic and non-hydrogenic QDT momentum space wavefunctions
;
; Paul Barklem 
; written Jan-March 2010
;
; notes:  1. QDT functions for l < 14 can be calculated easily with analytic functions,
;            but higher l uses series expansions which are very slow to compute
;            - if l > 13 needed, best to derive those Hankel transforms analytically
;         2. There may be phase differences of factor -1 between different calculations
;             - irrelevant for me since I want the square
;         3. Note limited validity of Coulomb approximation for gwf
;             - gwfsq derives the square of the wavefunction for all momenta by taking QDT
;               where valid and interpolating in hydrogenic wfs outside

function Fnl, A, ION, n, l, p
; calculates the hydrogenic wavefunction in momentum space 
; follows Bransden and Joachain, Physics of Atoms and Molecules, 2nd edition appendix 5
;
; where:
; A is the atomic mass
; ION is the total charge on the atom (0 = neutral, 1=singly ionised)
; n, l are usual quantum numbers
; p is momentum in au

ZZ = ION + 1                             ; charge on core (i.e atom minus Rydberg electron)
mu = A * 1836.15 / (A * 1836.15 + 1.)    ; reduced mass of the core-electron system in au
pmu = mu                                 ; p0 = 1 in au, multiply by mu/m to correct for reduced mass
ps = p/(ZZ*pmu)                          ; corrections for Z and mu
nd = double(n)
ld = double(l)
n2ps2 = nd*nd*ps*ps

part1 = sqrt(2.d0*factorial(n-l-1)/!pi/factorial(n+l))
part2 = nd*nd*2.d0^(2*ld+2)*factorial(l)
part3 = nd^ld*ps^ld/((n2ps2+1.d0)^(ld+2))
x = (n2ps2-1.d0)/(n2ps2+1.d0)
part4 = Gegenbauer(x,l+1,n-l-1)

return, part1 * part2 * part3 * part4
end

function gegenbauer, x, alpha, n
; calculates the Gegenbauer polynomial

if n lt 0 or n ne fix(n) then begin
   print, ' Gegenbauer: n<0 or non-integer not permitted'
   stop
endif

case n of
   0: C = 1.d0
   1: C = 2.d0 * alpha * x
   else: begin
       Cm2 = 1.d0
       Cm1 = 2.d0 * alpha * x
       for m = 2, n do begin
          C = (2.d0*x*(m+alpha-1.d0)*Cm1/m - (m+2.d0*alpha-2.d0)*Cm2/m) 
          Cm2 = Cm1
          Cm1 = C
       endfor
       end
endcase       

return, C
end

function hankel, q, nu, l, t
; the Hankel tranform we need by making appropriate substitutions to get
; eq 31 in Hoang Binh and van Regemorter
;
; for l < 14 we use analytic expressions calculated with Mathematica

case l of
   0 : x = sqrt(2./!pi) * gamma(nu-t+1) * ((1./nu/nu + q*q)^((-nu+t-1)/2.))*sin((nu-t+1)*atan(q*nu))
   1 : x = sqrt(2./!pi)/q * gamma(nu-t) * ((1./nu/nu + q*q)^((-nu+t)/2.))*sin((nu-t)*atan(q*nu)) $
         -sqrt(2./!pi) * gamma(nu-t+1) * ((1./nu/nu + q*q)^((-nu+t-1)/2.))*cos((nu-t+1)*atan(q*nu))
   2 : x =          -((nu^(-1 + nu - t)*Sqrt(2/!Pi)* $
            (1 + nu^2*q^2)^((-1 - nu + t)/2.)* $
            Gamma(-1 + nu - t)* $
            (3*nu*q*Sqrt(1 + nu^2*q^2)*(-1 + nu - t)* $
               Cos((-nu + t)*atan(nu*q)) + $
              nu^2*q^2*((-1 + nu)*nu + t - 2*nu*t + t^2)* $
               Sin((1 + nu - t)*atan(nu*q)) + $
              3*(1 + nu^2*q^2)*Sin((1 - nu + t)*atan(nu*q))))/ $
          q^2)
   3 :  x =         (nu^(-2 + nu - t)*Sqrt(2/!Pi)* $
          (1 + nu^2*q^2)^((-1 - nu + t)/2.)* $
          (-15*nu*q*(1 + nu^2*q^2)* $
             Cos((1 - nu + t)*atan(nu*q))*Gamma(-1 + nu - t) + $
            nu^3*q^3*Cos((1 + nu - t)*atan(nu*q))* $
             Gamma(1 + nu - t) + $
            3*Sqrt(1 + nu^2*q^2)* $
             (2*nu^2*q^2*Gamma(nu - t)* $
                Sin((-nu + t)*atan(nu*q)) - $
               5*(1 + nu^2*q^2)*Gamma(-2 + nu - t)* $
                Sin((2 - nu + t)*atan(nu*q)))))/q^3 
   4 : x =         (nu^(-3 + nu - t)*Sqrt(2/!Pi)* $
          (1 + nu^2*q^2)^((-1 - nu + t)/2.)* $
          (-105*nu*q*(1 + nu^2*q^2)^1.5* $
             Cos((2 - nu + t)*atan(nu*q))*Gamma(-2 + nu - t) + $
            10*nu^3*q^3*Sqrt(1 + nu^2*q^2)* $
             Cos((-nu + t)*atan(nu*q))*Gamma(nu - t) + $
            nu^4*q^4*Gamma(1 + nu - t)* $
             Sin((1 + nu - t)*atan(nu*q)) + $ 
            45*nu^2*q^2*Gamma(-1 + nu - t)* $
             Sin((1 - nu + t)*atan(nu*q)) + $
            45*nu^4*q^4*Gamma(-1 + nu - t)* $
             Sin((1 - nu + t)*atan(nu*q)) - $
            105*Gamma(-3 + nu - t)* $
             Sin((3 - nu + t)*atan(nu*q)) - $ 
            210*nu^2*q^2*Gamma(-3 + nu - t)* $
             Sin((3 - nu + t)*atan(nu*q)) - $
            105*nu^4*q^4*Gamma(-3 + nu - t)* $
             Sin((3 - nu + t)*atan(nu*q))))/q^4 
    5 : x =         -((nu^(-4 + nu - t)*Sqrt(2/!Pi)* $
            (1 + nu^2*q^2)^((-1 - nu + t)/2.)* $
            (945*nu*q*(1 + nu^2*q^2)^2* $
               Cos((3 - nu + t)*atan(nu*q))*Gamma(-3 + nu - t) - $
              105*nu^3*q^3*(1 + nu^2*q^2)* $
               Cos((1 - nu + t)*atan(nu*q))*Gamma(-1 + nu - t) + $
              nu^5*q^5*Cos((1 + nu - t)*atan(nu*q))* $
               Gamma(1 + nu - t) + $
              15*nu^4*q^4*Sqrt(1 + nu^2*q^2)*Gamma(nu - t)* $
               Sin((-nu + t)*atan(nu*q)) - $
              420*nu^2*q^2*Sqrt(1 + nu^2*q^2)* $
               Gamma(-2 + nu - t)*Sin((2 - nu + t)*atan(nu*q)) - $ 
              420*nu^4*q^4*Sqrt(1 + nu^2*q^2)* $
               Gamma(-2 + nu - t)*Sin((2 - nu + t)*atan(nu*q)) + $ 
              945*Sqrt(1 + nu^2*q^2)*Gamma(-4 + nu - t)* $
               Sin((4 - nu + t)*atan(nu*q)) + $
              1890*nu^2*q^2*Sqrt(1 + nu^2*q^2)* $
               Gamma(-4 + nu - t)*Sin((4 - nu + t)*atan(nu*q)) + $ 
              945*nu^4*q^4*Sqrt(1 + nu^2*q^2)* $
               Gamma(-4 + nu - t)*Sin((4 - nu + t)*atan(nu*q))))/ $
          q^5) 
    6 : x =         -((nu^(nu - t)*Sqrt(2/!Pi)* $
            (1 + nu^2*q^2)^((-nu + t)/2.)*Gamma(-5 + nu - t)* $
            (-21*q^5*(1 - nu + t)*(2 - nu + t)*(3 - nu + t)* $
               (4 - nu + t)*(5 - nu + t)* $
               Cos((-nu + t)*atan(nu*q)) + $ 
              (1260*q^3*(1 + nu^2*q^2)*(3 - nu + t)* $
                 (4 - nu + t)*(5 - nu + t)* $
                 Cos((2 - nu + t)*atan(nu*q)))/nu^2 + $ 
              (10395*q*(1 + nu^2*q^2)^2*(-5 + nu - t)* $
                 Cos((4 - nu + t)*atan(nu*q)))/nu^4 + $
              (nu*q^6*(-5 + nu - t)*(-4 + nu - t)*(-3 + nu - t)* $
                 (-2 + nu - t)*(-1 + nu - t)*(nu - t)* $
                 Sin((1 + nu - t)*atan(nu*q)))/ $
               Sqrt(1 + nu^2*q^2) + $
              (210*q^4*Sqrt(1 + nu^2*q^2)*(-5 + nu - t)* $
                 (-4 + nu - t)*(-3 + nu - t)*(-2 + nu - t)* $
                 Sin((1 - nu + t)*atan(nu*q)))/nu - $
              (4725*q^2*(1 + nu^2*q^2)^1.5*(-5 + nu - t)* $
                 (-4 + nu - t)*Sin((3 - nu + t)*atan(nu*q)))/ $
               nu^3 + (10395*(1 + nu^2*q^2)^2.5* $
                 Sin((5 - nu + t)*atan(nu*q)))/nu^5))/q^6)
    7 : x =         (nu^(nu - t)*Sqrt(2/!Pi)*(1 + nu^2*q^2)^((-nu + t)/2.)*$
          ((-135135*q*(1 + nu^2*q^2)^2.5*$
               Cos((5 - nu + t)*atan(nu*q))*Gamma(-5 + nu - t))/$
             nu^5 + (17325*q^3*(1 + nu^2*q^2)^1.5*$
               Cos((3 - nu + t)*atan(nu*q))*Gamma(-3 + nu - t))/$
             nu^3 - (378*q^5*Sqrt(1 + nu^2*q^2)*$
               Cos((1 - nu + t)*atan(nu*q))*Gamma(-1 + nu - t))/$
             nu + (nu*q^7*Cos((1 + nu - t)*atan(nu*q))*$
               Gamma(1 + nu - t))/Sqrt(1 + nu^2*q^2) + $
            28*q^6*Gamma(nu - t)*Sin((-nu + t)*atan(nu*q)) - $
            (3150*q^4*(1 + nu^2*q^2)*Gamma(-2 + nu - t)*$
               Sin((2 - nu + t)*atan(nu*q)))/nu^2 + $
            (62370*(q + nu^2*q^3)^2*Gamma(-4 + nu - t)*$
               Sin((4 - nu + t)*atan(nu*q)))/nu^4 - $
            (135135*(1 + nu^2*q^2)^3*Gamma(-6 + nu - t)*$
               Sin((6 - nu + t)*atan(nu*q)))/nu^6))/q^7 
    8 : x =         (nu^(nu - t)*Sqrt(2/!Pi)*(1 + nu^2*q^2)^((-nu + t)/2.)*$
          ((-2027025*q*(1 + nu^2*q^2)^3*$
               Cos((6 - nu + t)*atan(nu*q))*Gamma(-6 + nu - t))/$
             nu^6 + (270270*q^3*(1 + nu^2*q^2)^2*$
               Cos((4 - nu + t)*atan(nu*q))*Gamma(-4 + nu - t))/$
             nu^4 - (6930*q^5*(1 + nu^2*q^2)*$
               Cos((2 - nu + t)*atan(nu*q))*Gamma(-2 + nu - t))/$
             nu^2 + 36*q^7*Cos((-nu + t)*atan(nu*q))*$
             Gamma(nu - t) + $
            (nu*q^8*Gamma(1 + nu - t)*$
               Sin((1 + nu - t)*atan(nu*q)))/Sqrt(1 + nu^2*q^2)$
              + (630*q^6*Sqrt(1 + nu^2*q^2)*Gamma(-1 + nu - t)*$
               Sin((1 - nu + t)*atan(nu*q)))/nu - $
            (51975*q^4*(1 + nu^2*q^2)^1.5*Gamma(-3 + nu - t)*$
               Sin((3 - nu + t)*atan(nu*q)))/nu^3 + $
            (945945*q^2*(1 + nu^2*q^2)^2.5*Gamma(-5 + nu - t)*$
               Sin((5 - nu + t)*atan(nu*q)))/nu^5 - $
            (2027025*(1 + nu^2*q^2)^3.5*Gamma(-7 + nu - t)*$
               Sin((7 - nu + t)*atan(nu*q)))/nu^7))/q^8
    9 : x =         (nu^(nu - t)*Sqrt(2/!Pi)*(1 + nu^2*q^2)^((-nu + t)/2.)*$
          ((-34459425*q*(1 + nu^2*q^2)^3.5*$
               Cos((7 - nu + t)*atan(nu*q))*Gamma(-7 + nu - t))/$
             nu^7 + (4729725*q^3*(1 + nu^2*q^2)^2.5*$
               Cos((5 - nu + t)*atan(nu*q))*Gamma(-5 + nu - t))/$
             nu^5 - (135135*q^5*(1 + nu^2*q^2)^1.5*$
               Cos((3 - nu + t)*atan(nu*q))*Gamma(-3 + nu - t))/$
             nu^3 + (990*q^7*Sqrt(1 + nu^2*q^2)*$
               Cos((1 - nu + t)*atan(nu*q))*Gamma(-1 + nu - t))/$
             nu - (nu*q^9*Cos((1 + nu - t)*atan(nu*q))*$
               Gamma(1 + nu - t))/Sqrt(1 + nu^2*q^2) - $
            45*q^8*Gamma(nu - t)*Sin((-nu + t)*atan(nu*q)) + $
            (13860*q^6*(1 + nu^2*q^2)*Gamma(-2 + nu - t)*$
               Sin((2 - nu + t)*atan(nu*q)))/nu^2 - $
            (945945*q^4*(1 + nu^2*q^2)^2*Gamma(-4 + nu - t)*$
               Sin((4 - nu + t)*atan(nu*q)))/nu^4 + $
            (16216200*q^2*(1 + nu^2*q^2)^3*Gamma(-6 + nu - t)*$
               Sin((6 - nu + t)*atan(nu*q)))/nu^6 - $
            (34459425*(1 + nu^2*q^2)^4*Gamma(-8 + nu - t)*$
               Sin((8 - nu + t)*atan(nu*q)))/nu^8))/q^9
    10 :x =         -((nu^(nu - t)*Sqrt(2/!Pi)* $
            (1 + nu^2*q^2)^((-nu + t)/2.)*$
            ((654729075*q*(1 + nu^2*q^2)^4*$
                 Cos((8 - nu + t)*atan(nu*q))*Gamma(-8 + nu - t))$
                /nu^8 - (91891800*(q + nu^2*q^3)^3*$
                 Cos((6 - nu + t)*atan(nu*q))*Gamma(-6 + nu - t))$
                /nu^6 + (2837835*q^5*(1 + nu^2*q^2)^2*$
                 Cos((4 - nu + t)*atan(nu*q))*Gamma(-4 + nu - t))$
                /nu^4 - (25740*q^7*(1 + nu^2*q^2)*$
                 Cos((2 - nu + t)*atan(nu*q))*Gamma(-2 + nu - t))$
                /nu^2 + 55*q^9*Cos((-nu + t)*atan(nu*q))*$
               Gamma(nu - t) + $
              (nu*q^10*Gamma(1 + nu - t)*$
                 Sin((1 + nu - t)*atan(nu*q)))/$
               Sqrt(1 + nu^2*q^2) + $
              (1485*q^8*Sqrt(1 + nu^2*q^2)*Gamma(-1 + nu - t)*$
                 Sin((1 - nu + t)*atan(nu*q)))/nu - $
              (315315*q^6*(1 + nu^2*q^2)^1.5*$
                 Gamma(-3 + nu - t)*Sin((3 - nu + t)*atan(nu*q)))$
                /nu^3 + (18918900*q^4*(1 + nu^2*q^2)^2.5*$
                 Gamma(-5 + nu - t)*Sin((5 - nu + t)*atan(nu*q)))$
                /nu^5 - (310134825*q^2*(1 + nu^2*q^2)^3.5*$
                 Gamma(-7 + nu - t)*Sin((7 - nu + t)*atan(nu*q)))$
                /nu^7 + (654729075*(1 + nu^2*q^2)^4.5*$
                 Gamma(-9 + nu - t)*Sin((9 - nu + t)*atan(nu*q)))$
               /nu^9))/q^10) 
    11 :x =         (nu^(nu - t)*Sqrt(2/!Pi)*(1 + nu^2*q^2)^((-nu + t)/2.)*$
          ((-13749310575*q*(1 + nu^2*q^2)^4.5*$
               Cos((9 - nu + t)*atan(nu*q))*Gamma(-9 + nu - t))/$
             nu^9 + (1964187225*q^3*(1 + nu^2*q^2)^3.5*$
               Cos((7 - nu + t)*atan(nu*q))*Gamma(-7 + nu - t))/$
             nu^7 - (64324260*q^5*(1 + nu^2*q^2)^2.5*$
               Cos((5 - nu + t)*atan(nu*q))*Gamma(-5 + nu - t))/$
             nu^5 + (675675*q^7*(1 + nu^2*q^2)^1.5*$
               Cos((3 - nu + t)*atan(nu*q))*Gamma(-3 + nu - t))/$
             nu^3 - (2145*q^9*Sqrt(1 + nu^2*q^2)*$
               Cos((1 - nu + t)*atan(nu*q))*Gamma(-1 + nu - t))/$
             nu + (nu*q^11*Cos((1 + nu - t)*atan(nu*q))*$
               Gamma(1 + nu - t))/Sqrt(1 + nu^2*q^2) + $
            66*q^10*Gamma(nu - t)*Sin((-nu + t)*atan(nu*q)) - $
            (45045*q^8*(1 + nu^2*q^2)*Gamma(-2 + nu - t)*$
               Sin((2 - nu + t)*atan(nu*q)))/nu^2 + $
            (7567560*q^6*(1 + nu^2*q^2)^2*Gamma(-4 + nu - t)*$
               Sin((4 - nu + t)*atan(nu*q)))/nu^4 - $
            (413513100*q^4*(1 + nu^2*q^2)^3*Gamma(-6 + nu - t)*$
               Sin((6 - nu + t)*atan(nu*q)))/nu^6 + $
            (6547290750*q^2*(1 + nu^2*q^2)^4*$
               Gamma(-8 + nu - t)*Sin((8 - nu + t)*atan(nu*q)))/$
             nu^8 - (13749310575*(1 + nu^2*q^2)^5*$
               Gamma(-10 + nu - t)*Sin((10 - nu + t)*atan(nu*q)))$
              /nu^10))/q^11         
    12 :x =         (Sqrt(2/!Pi)*(-316234143225*(1/nu)^(-nu + t)*q*$
             (nu^(-2) + q^2)^5*(1 + nu^2*q^2)^((-nu + t)/2.)*$
             Cos(10*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-10 + nu - t) + $
            45831035250*(1/nu)^(-nu + t)*q^3*$
              (nu^(-2) + q^2)^4*(1 + nu^2*q^2)^((-nu + t)/2.)*$
             Cos(8*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-8 + nu - t) - $
            1571349780*(1/nu)^(-nu + t)*q^5*(nu^(-2) + q^2)^3*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*$
             Cos(6*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-6 + nu - t) + $
            18378360*(1/nu)^(-nu + t)*q^7*(nu^(-2) + q^2)^2*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*$
             Cos(4*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-4 + nu - t) - $
            75075*(1/nu)^(-nu + t)*q^9*(nu^(-2) + q^2)*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*$
             Cos(2*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-2 + nu - t) + $
            78*(1/nu)^(-nu + t)*q^11*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*$
             Cos(nu*atan(nu*q) - t*atan(nu*q))*Gamma(nu - t) + $
            (1/nu)^(-1 - nu + t)*q^12*$
             (1 + nu^2*q^2)^((-1 - nu + t)/2.)*$
             Gamma(1 + nu - t)*$
             Sin(atan(nu*q) + nu*atan(nu*q) - t*atan(nu*q)) $
             + 3003*(1/nu)^(1 - nu + t)*q^10*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Gamma(-1 + nu - t)*$
             Sin(atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q)) $
             - 1351350*(1/nu)^(1 - nu + t)*q^8*(nu^(-2) + q^2)*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Gamma(-3 + nu - t)*$
             Sin(3*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              + 192972780*(1/nu)^(1 - nu + t)*q^6*$
             (nu^(-2) + q^2)^2*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Gamma(-5 + nu - t)*$
             Sin(5*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              - 9820936125*(1/nu)^(1 - nu + t)*q^4*$
             (nu^(-2) + q^2)^3*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Gamma(-7 + nu - t)*$
             Sin(7*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              + 151242416325*(1/nu)^(1 - nu + t)*q^2*$
             (nu^(-2) + q^2)^4*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Gamma(-9 + nu - t)*$
             Sin(9*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              - 316234143225*(1/nu)^(1 - nu + t)*$
             (nu^(-2) + q^2)^5*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Gamma(-11 + nu - t)*$
             Sin(11*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))))/q^12
    13: x =         (Sqrt(2/!Pi)*(-7905853580625*(1/nu)^(1 - nu + t)*q*$
             (nu^(-2) + q^2)^5*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Cos(11*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-11 + nu - t) + $
            1159525191825*(1/nu)^(1 - nu + t)*q^3*$
             (nu^(-2) + q^2)^4*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Cos(9*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-9 + nu - t) - $
            41247931725*(1/nu)^(1 - nu + t)*q^5*$
             (nu^(-2) + q^2)^3*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Cos(7*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-7 + nu - t) + $
            523783260*(1/nu)^(1 - nu + t)*q^7*$
             (nu^(-2) + q^2)^2*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Cos(5*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-5 + nu - t) - $
            2552550*(1/nu)^(1 - nu + t)*q^9*(nu^(-2) + q^2)*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Cos(3*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q))*Gamma(-3 + nu - t) + $
            4095*(1/nu)^(1 - nu + t)*q^11*$
             (1 + nu^2*q^2)^((1 - nu + t)/2.)*$
             Cos(atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))*$
             Gamma(-1 + nu - t) - $
            (1/nu)^(-1 - nu + t)*q^13*$
             (1 + nu^2*q^2)^((-1 - nu + t)/2.)*$
             Cos(atan(nu*q) + nu*atan(nu*q) - t*atan(nu*q))*$
             Gamma(1 + nu - t) + $
            91*(1/nu)^(-nu + t)*q^12*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*Gamma(nu - t)*$
             Sin(nu*atan(nu*q) - t*atan(nu*q)) + $
            120120*(1/nu)^(-nu + t)*q^10*(nu^(-2) + q^2)*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*Gamma(-2 + nu - t)*$
             Sin(2*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              - 41351310*(1/nu)^(-nu + t)*q^8*$
             (nu^(-2) + q^2)^2*(1 + nu^2*q^2)^((-nu + t)/2.)*$
             Gamma(-4 + nu - t)*$
             Sin(4*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              + 5237832600*(1/nu)^(-nu + t)*q^6*$
             (nu^(-2) + q^2)^3*(1 + nu^2*q^2)^((-nu + t)/2.)*$
             Gamma(-6 + nu - t)*$
             Sin(6*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              - 252070693875*(1/nu)^(-nu + t)*q^4*$
             (nu^(-2) + q^2)^4*(1 + nu^2*q^2)^((-nu + t)/2.)*$
             Gamma(-8 + nu - t)*$
             Sin(8*atan(nu*q) - nu*atan(nu*q) + t*atan(nu*q))$
              + 3794809718700*(1/nu)^(-nu + t)*q^2*$
             (nu^(-2) + q^2)^5*(1 + nu^2*q^2)^((-nu + t)/2.)*$
             Gamma(-10 + nu - t)*$
             Sin(10*atan(nu*q) - nu*atan(nu*q) + $
               t*atan(nu*q)) - $
            7905853580625*(1/nu)^(-nu + t)*(nu^(-2) + q^2)^6*$
             (1 + nu^2*q^2)^((-nu + t)/2.)*Gamma(-12 + nu - t)*$
             Sin(12*atan(nu*q) - nu*atan(nu*q) + $
              t*atan(nu*q))))/q^13
    else: begin
          print, 'llow > 13: slow calculation - not recommended!'
          a = 1./double(nu)
          b = double(q)
          nup = l+0.5
          mu = double(nu) - t + 1.5
          x = watson(a,b,nup,mu) * sqrt(b) 
         end
endcase
return, x
end

function hyper2F1, a, b, c, z
; calculates the 2F1 hypergeometric series

if (a eq 0.) or (b eq 0.) or (c eq 0.) then return, 0.d0

errlim = 1d-16    
err = 1d10
sum = 0.d0
n = 0l
np = fix(6000*z)>100 ; number of terms per iteration, some tests showed this quite good 
           
prodlast = 1.d0              
while (err gt errlim) do begin
  nn = n + indgen(np)
  prod = nn * 0.d0 + 1.d0
  iarr = [0.d0, dindgen(nn[np-1])]
  iarr = iarr[nn[0]:nn[np-1]]
  ; term by term evaluation of product for numerical precision
  ; product is accumulated for speed
  parr = (a+iarr) / (c+iarr) * (b+iarr) / (iarr+1.) 
  prod[0] = prodlast
  if nn[0] ge 1 then prod[0] =  prodlast*parr[0]
  for i = 1, np-1 do prod[i] =  prod[i-1]*parr[i]   
  term = prod * double(z)^nn
  if n gt 1000 then err = abs(total(term))/abs(sum) else err = min(abs(term))/abs(sum)
  sum = sum + total(term)
  n = n + np
  prodlast = prod[np-1]
endwhile

return, sum
end

function watson, a, b, nu, mu
; general Hankel transform int(exp(-at) J_nu(bt) t^(mu-1) dt, 0, infinity)
; pg 385, Watson - Theory of Bessel functions
;
; there are alternate expressions for this given by Watson
; this expression is chosen since z < 1 always.  
; This means the 2F1 should converge 

ad = double(a)
bd = double(b)
nud= double(nu)
mud = double(mu)

asqbsq = ad*ad+bd*bd
nuplmu = nud + mud

term1 = (0.5*bd)^nud * gamma(nuplmu) / (asqbsq)^(0.5*(nuplmu)) / gamma(nud+1)
term2 = hyper2F1((nuplmu)/2., (1.-mud+nud)/2., nud+1, bd*bd/(asqbsq))

return, term1*term2
end


function Knorm, nu, l
; normalisation factor from Hoang Binh and van Regemorter, eq 15
; assuming eta = 1, true if nu is not small

nud = double(nu)
return, 1.d0/(nud*sqrt(gamma(nud+l+1)*gamma(nud-l)))
end

function bfunc, nu, l, t
; calculate the b coefficients
; the expressions of Bates and Damgaard 1949 are used since they work everywhere (even at integer nu values)
; if one doesn't have small, integer nu values then one could use:
; b = (-1.d0)^td * gamma(-nu+l+1+td)*gamma(-nu-l+td)/gamma(-nu+l+1)/gamma(-nu-l)/factorial(td) *(nu/2.d0)^td 
; from Hoang Binh and van Regemorter, but note (nu/2)^t term is missing from paper
;
; The two methods seem to be equivalent to very good precision even at large nu and l

td = double(t)
if t eq 0 then begin
   b = 1.
endif
if t ge 1 then begin
   blast = 1.
   for i = 1, t do begin
   b = -nu/2./i * (nu-l-i)*(nu+l+1.-i) * blast
   blast = b
   endfor
endif

return, b
end

function gwf, nu, l, q
; the QDT non-hydrogenic wavefunction in momentum space, following Hoang Binh and van Regemorter (1997)
; where: 
; nu is the effective principal quantum number
; l  is angular momentum quantum number
; q  is momentum in atomic units
; normalisation is Integral |g|^2 q^2 dq = 1
;
; note due to machine numerical precision limits, if nu > 24, we use the nearest hydrogenic wavefunction

nud = double(nu)
nq = n_elements(q)
wf = dblarr(nq)

if (nu gt 20.) or (abs(nu mod 1) lt 0.001) then begin     ; numerical precision limits QDT wf calc to about n*~25
                                                          ; also use exact analytic expression if integer nu value (n*)
   wf = Fnl(1., 0, round(nu), l, q)
endif else begin

t0 = fix(nud + l + 1.) 

for i = 0, nq-1 do begin
   if q[i] le 0.d0 then begin
      wf[i] = 0.d0  
   endif else begin
      term1 = Knorm(nud, l)/q[i] * (2.d0/nud)^nud
      sum = 0.d0
      for t = 0, t0 do begin
         term2 = bfunc(nud, l, t) * hankel(q[i], nud, l, t) 
         sum = sum + term2 
      endfor
      wf[i] = sum * term1 
   endelse
endfor
endelse

return, wf
end

function gwfsq, nu, l, q
; generate a best estimate of square of wavefunction 

nud = double(nu)
nq = n_elements(q)
wf = dblarr(nq)

; set bound for using Coulomb approx QDT functions based on empirical study accounting for
; breakdown of Coulomb approx, loss of numerical precision, and where wf has smooth behaviour
; and can be interpolated
; these cutoffs have been tested for 2.5 < nu < 20.5 and l <= 13

qbar = 2/!pi/nud
qmin = qbar
if nu le 11. then qmin = qbar / 2.
if nu le 8. then qmin = qbar / 4.
qmax = (nu - l + 1) * qbar  ; this seems to work ok.
if nu ge 8.5 then qmax = ((nu - l)>2) * qbar
if nu ge 11.5 then qmax = ((nu - l - 1)/2.>2) * qbar

nup = fix(nud+1)   ; nearest integers n
num = fix(nud)
if num eq l then begin
   num = num+1
   nup = nup+1
endif

for i = 0, nq-1 do begin
   if q[i] le 0.d0 then begin
      wf[i] = 0.d0 
   endif else begin
      if (q[i] lt qmax) and (q[i] gt qmin) then begin  ; use QDT
         wf[i] = gwf(nud, l, q[i])^2.
      endif else begin                                 ; linear interpolation between integer values for nu
         ym = Fnl(1., 0, num, l, q[i])^2.
         yp = Fnl(1., 0, nup, l, q[i])^2.
         h = (nu - num)/(nup - num)
         wf[i] = 10^(alog10(ym) + h*(alog10(yp) - alog10(ym)))
      endelse
   endelse
endfor

return, wf
end




