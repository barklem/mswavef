
function hyper2F1old, a, b, c, z
; calculates the 2F1 hypergeometric series
; this is an older version, non-vectorised version, but well tested.

if (a eq 0.) or (b eq 0.) or (c eq 0.) then return, 0.d0

errlim = 1d-16  
err = 1d10
sum = 0.d0
n = 0l
while (err gt errlim) do begin
  prod = 1.d0 
  ; term by term evaluation of product for numerical precision  
  if n ge 1 then begin
     iarr = dindgen(n)
     prod =  product(1.d0 / (iarr+1.) * (a+iarr) / (c+iarr) * (b+iarr))
  endif
  term = prod * double(z)^n
  if n ge 1 then err = abs(term/sum)
  sum = sum + term
  n = n + 1
endwhile

return, sum
end


pro testhankel
q = dindgen(100)/100.+1
nu = 10.0
l = 9
t = 4
!p.multi=[0,1,3]
y1 = q*0.d0
for i = 0, n_elements(q)-1 do y1[i] = hankel(q[i],nu,l,t)
y2 = sqrt(2./!pi) * gamma(nu-t+1) * ((1./nu/nu + q*q)^((-nu+t-1)/2.))*sin((nu-t+1)*atan(q*nu))
plot, q, y1, charsize = 2
plot, q, y2, psym = 5, charsize = 2
plot, q, y1/y2, charsize = 2, thicK = 3
!p.multi=0
stop
end

pro testhyper

z = dindgen(100)/100.
a = .5
b = 1.
c = 1.5
h1 = z*0.d0
for i = 0, n_elements(z)-1 do h1[i] = z[i] * hyper2F1(a,b,c,z[i]*z[i])
h2 = 0.5*alog((1+z)/(1-z))
plot, z, h2, linestyle=2, thick = 3
oplot, z, h1
stop
end

pro testgwf
t = systime(1)
q = (dindgen(100))/40
window, 0, xsize = 1900
!p.multi = [0,5,1]
!p.charsize = 3
; Na 4s
ns = 2.64
plot, q, q*gwf(ns, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)', title = 'Na (4s)'
oplot, q, q*gwf(2.0, 0, q), linestyle = 2
oplot, q, q*gwf(3.0, 0, q), linestyle = 2
;oplot, q, (q*gwf(2.2, 0, q))^2, linestyle = 3
;oplot, q, (q*gwf(2.8, 0, q))^2, linestyle = 3
oplot, q, q*Fnl(1., 0, 3, 0, q), psym = 4
; Na 4p
ns = 3.13
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = 'Na (4p)'
oplot, q, q*gwf(3.0, 1, q), linestyle = 2
oplot, q, -q*Fnl(1., 0, 3, 1, q), psym = 4
; Mg 6s
ns = 4.358
plot, q, q*gwf(ns, 0, q), title = 'Mg (6s)'
oplot, q, q*gwf(4.0, 2, q), linestyle = 2
oplot, q, -q*Fnl(1., 0., 4.0, 2, q), psym = 4
; Mg 6p
ns = 4.86
plot, q, q*gwf(ns, 1, q), title = 'Mg (6p)';, xr=[0,1.0], /xs
oplot, q, q*gwf(5.0, 1, q), linestyle = 2
oplot, q, -q*Fnl(1., 0, 5, 1, q), psym = 4
; Cs 10p
ns = 6.43
plot, q, q*gwf(ns, 1, q), title = 'Cs (10p)';, xr=[0,1.0], /xs
oplot, q, q*gwf(6.0, 1, q), linestyle = 2
oplot, q, q*Fnl(1., 0, 6, 1, q), psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwfl
t = systime(1)
q = [(dindgen(20)+1), 50.]
;window, 0, xsize = 1700, ysize = 1000
set_plot, 'ps'
device, file = 'wfs.eps', ysize = 30, xsize = 45
!p.multi = [0,6,4]
!p.charsize = 2
; Na 4s
l = 0
ns = 2.1
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 5.7
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs 
print, systime(1) - t, ' SECONDS'
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'

l = 2
ns = 5.7
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 10.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS' 
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'


l = 4
ns = 5.7
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 10.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS' 
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'

l = 7
ns = 8.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 10.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS' 
ns = 12.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
print, systime(1) - t, ' SECONDS'

l = 13
ns = 14.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
ns = 15.0
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs 
ns = 16.5
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
ns = 17.3
yy = gwf(ns, l, q)
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (q*yy)^2, ytitle = "q*g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /xlog, xr=[min(q),max(q)], /xs
plot, q, (yy)^2, ytitle = "g(q)", xtitle = 'q (au)', /ylog, xr=[min(q),max(q)], /xs
;plot, q, (q*gwf(2.2, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*gwf(2.8, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*gwf(2.0, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*gwf(3.0, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
;oplot, q, (q*Fnl(1., 0, 3, 0, q))^2, psym = 4
print, systime(1) - t, ' SECONDS'
device, /close_file
set_plot, 'x'
stop
; Na 4p
ns = 3.13
plot, q, (gwf(ns, 1, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
oplot, q, (gwf(3.0, 1, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 3, 1, q))^2, psym = 4
; Mg 6s
ns = 4.358
plot, q, (gwf(ns, 0, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs
oplot, q, (gwf(4.0, 0, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 4, 0, q))^2, psym = 4
; Mg 6p
ns = 4.86
plot, q, (gwf(ns, 1, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs ;, xr=[0,1.0], /xs
oplot, q, (gwf(5.0, 1, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 5, 1, q))^2, psym = 4
; Cs 10p
ns = 6.43
plot, q, (gwf(ns, 1, q))^2, /ylog, /xlog, xr=[min(q),max(q)], /xs ;, xr=[0,1.0], /xs
oplot, q, (gwf(6.0, 1, q))^2, linestyle = 2
oplot, q, (Fnl(1., 0, 6, 1, q))^2, psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwfl2
t = systime(1)
pmin = 1e-4
pmax = 1000.
npts = 10000
pstep = (alog(pmax) - alog(pmin))/(npts-1)
q = exp(alog(pmin) + pstep*dindgen(npts))
set_plot, 'ps'
device, file = 'wfs.eps', xsize = 28, xoffset=1, yoffset = 29, /landscape
!p.multi = 0
!p.charsize = 1.5

ns = 20.0
yy = gwf(ns, 0, q)
plot, q, (yy)^2, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')
for ns = 3.0, 20.0, 1.0 do begin
yy = gwf(ns, 0, q)
oplot, q, (yy)^2, linestyle = 5
endfor

for ns = 3.0, 20.0, 1.0 do begin

yy = gwf(ns, 0, q)
plot, q, (yy)^2, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')

for l = 1, min([13, fix(ns)-1]), 1 do begin
oplot, q, gwf(ns, l, q)^2, linestyle = 5
end

for l = 0, min([13, fix(ns)-1]), 1 do begin

yy = gwf(ns, l, q)
plot, q, (yy)^2, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
oplot, q, gwf(ns-0.5, l, q)^2, linestyle = 5
oplot, q, gwf(ns-0.4, l, q)^2,linestyle = 5
oplot, q, gwf(ns-0.3, l, q)^2, linestyle = 5
oplot, q, gwf(ns-0.2, l, q)^2, linestyle = 5
oplot, q, gwf(ns-0.1, l, q)^2, linestyle = 5
oplot, q, gwf(ns+0.5, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.4, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.3, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.2, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.1, l, q)^2, linestyle = 1
oplot, q, gwf(ns+0.01, l, q)^2, linestyle = 1, thick = 10

oplot, [2/!Pi/ns,2/!Pi/ns], [1e-30,1e30] 
oplot, [2/!Pi/ns,2/!Pi/ns]*2, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*4, [1e-30,1e30], linestyle = 1 
oplot, [2/!Pi/ns,2/!Pi/ns]*.5, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*.25, [1e-30,1e30], linestyle = 1 
print, systime(1) - t, ' SECONDS'
end
end

device, /close_file
set_plot, 'x'
!p.multi = 0
!p.charsize = 1
end

pro testgwfl3
t = systime(1)
pmin = 1e-4
pmax = 1000.
npts = 10000
pstep = (alog(pmax) - alog(pmin))/(npts-1)
q = exp(alog(pmin) + pstep*dindgen(npts))
set_plot, 'ps'
device, file = 'wfs_new.eps', xsize = 28, xoffset=1, yoffset = 29, /landscape
!p.multi = 0
!p.charsize = 1.5

ns = 20.0
yy = gwfsq(ns, 0, q)
plot, q, yy, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')
for ns = 3.0, 20.0, 1.0 do begin
yy = gwfsq(ns, 0, q)
oplot, q, yy, linestyle = 5
endfor

for ns = 3.0, 20.0, 1.0 do begin

yy = gwfsq(ns, 0, q)
plot, q, yy, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')

for l = 0, min([13, fix(ns)-1]), 1 do begin
oplot, q, gwfsq(ns, l, q), linestyle = 5
end

for l = 0, min([13, fix(ns)-1]), 1 do begin

print, ns, l

yy = gwfsq(ns, l, q)
plot, q, yy, ytitle = "g(q)^2", xtitle = 'q (au)', /ylog, /xlog, xr=[min(q),max(q)], /xs, title = string(ns, '(f4.1)')+string(l, '(i3)')
oplot, q, gwfsq(ns-0.5, l, q), linestyle = 5
oplot, q, gwfsq(ns-0.4, l, q),linestyle = 5
oplot, q, gwfsq(ns-0.3, l, q), linestyle = 5
oplot, q, gwfsq(ns-0.2, l, q), linestyle = 5
oplot, q, gwfsq(ns-0.1, l, q), linestyle = 5
oplot, q, gwfsq(ns+0.5, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.4, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.3, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.2, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.1, l, q), linestyle = 1
oplot, q, gwfsq(ns+0.01, l, q), linestyle = 1, thick = 10

oplot, [2/!Pi/ns,2/!Pi/ns], [1e-30,1e30] 
oplot, [2/!Pi/ns,2/!Pi/ns]*2, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*4, [1e-30,1e30], linestyle = 1 
oplot, [2/!Pi/ns,2/!Pi/ns]*.5, [1e-30,1e30], linestyle = 5 
oplot, [2/!Pi/ns,2/!Pi/ns]*.25, [1e-30,1e30], linestyle = 1 
print, systime(1) - t, ' SECONDS'
end
end

device, /close_file
set_plot, 'x'
!p.multi = 0
!p.charsize = 1
end


pro testgwf2
t = systime(1)
q = (dindgen(1000))/10000.
window, 0, xsize = 1900
!p.multi = [0,5,1]
!p.charsize = 3
ns = 20.
plot, q, q*gwf(ns, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)', title = '20s'
ns = 20.
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = '20p'
ns = 25.
plot, q, q*gwf(ns, 0, q), xtitle = 'q (au)', title = '25s'
ns = 25.
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = '25p';, xr=[0,1.0], /xs
ns = 30.
plot, q, q*gwf(ns, 1, q), xtitle = 'q (au)', title = '30p';, xr=[0,1.0], /xs
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwf3
t = systime(1)
print, systime()
q = (dindgen(10))/10000. + 0.0495001
q2 = (dindgen(1000))/1000000. + 0.0495
window, 0, xsize = 1900
!p.multi = [0,2,1]
!p.charsize = 3
ns = 10.
y1 = q*gwf(ns, 0, q)
y1a = -q2*Fnl(1., 0, 10, 0, q2)
plot, q, y1, ytitle = "q*g(q)", xtitle = 'q (au)', title = '20s'
oplot, q, y1a, psym = 4
print, systime(1) - t, ' SECONDS'
y2 = q*gwf(ns, 1, q)
y2a = q2*Fnl(1., 0, 20, 1, q2)
plot, q, y2, xtitle = 'q (au)', title = '20p'
oplot, q2, y2a, psym = 4
print, systime(1) - t, ' SECONDS'
print, systime()
!p.multi = 0
!p.charsize = 1
end


pro testgwf4
t = systime(1)
q = (dindgen(200))/1000.
window, 0, xsize = 1900
!p.multi = [0,5,1]
!p.charsize = 3
plot, q, q*gwf(15, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)'
oplot, q, q*Fnl(1., 0, 15, 0, q), psym = 4

plot, q, q*gwf(20, 0, q)
oplot, q, -q*Fnl(1., 0, 20, 0, q), psym = 4

plot, q, q*gwf(25, 0, q)
oplot, q, q*Fnl(1., 0, 25, 0, q), psym = 4

plot, q, q*gwf(30, 0, q)
oplot, q, -q*Fnl(1., 0, 30, 0, q), psym = 4

plot, q, q*gwf(50, 0, q)
oplot, q, q*Fnl(1., 0, 50, 0, q), psym = 4

print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwf5
t = systime(1)
q = (dindgen(1000))/2000.
window, 0, xsize = 1900
!p.multi = [0,1,1]
!p.charsize = 3
plot, q, q*gwf(9, 0, q), yr = [-4,4]
oplot, q, q*Fnl(1., 0, 9, 0, q), psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end

pro testgwf6
t = systime(1)
q = (dindgen(500))/5000.
window, 0, xsize = 1900
!p.multi = [0,1,1]
!p.charsize = 3
plot, q, q*Fnl(1., 0, 30, 0, q)
;oplot, q, q*Fnl(1., 0, 30, 0, q), linestyle = 1
oplot, q, -q*gwf(30., 0, q), linestyle = 5
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end
