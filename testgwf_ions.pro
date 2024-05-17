pro testgwf_ions
;
; testing for ions
;

t = systime(1)
q = (dindgen(200))/40
window, 0, xsize = 1900
!p.multi = [0,3,1]
!p.charsize = 3
; Na 4s
ns = 2.64
plot, q, q*gwf(0, ns, 0, q), ytitle = "q*g(q)", xtitle = 'q (au)', title = 'Na (4s)'
oplot, q, q*gwf(1, ns, 0, q), linestyle = 2
oplot, q, q*gwf(1, 3, 0, q), linestyle = 3
oplot, q, q*gwf(1, 2, 0, q), linestyle = 4
oplot, q, q*Fnl(1., 1, 3, 0, q), psym = 4

plot, q, (q*gwf(0, ns, 0, q))^2, ytitle = "(q*g(q))^2", xtitle = 'q (au)', title = 'Na (4s)'
oplot, q, (q*gwf(1, ns, 0, q))^2, linestyle = 2
oplot, q, (q*gwf(1, 3, 0, q))^2, linestyle = 3
oplot, q, (q*gwf(1, 2, 0, q))^2, linestyle = 4
oplot, q, (q*Fnl(1., 1, 3, 0, q))^2, psym = 4

plot, q, (gwf(0, ns, 0, q))^2, ytitle = "(g(q))^2", xtitle = 'q (au)', title = 'Na (4s)', /ylog
oplot, q, (gwf(1, ns, 0, q))^2, linestyle = 2
oplot, q, (gwf(1, 3, 0, q))^2, linestyle = 3
oplot, q, (gwf(1, 2, 0, q))^2, linestyle = 4
oplot, q, (Fnl(1., 1, 3, 0, q))^2, psym = 4


print, int_tabulated(q, (q*gwf(0, ns, 0, q))^2)
print, int_tabulated(q, (q*gwf(1, ns, 0, q))^2)

stop


; Na 4p
ns = 3.13
plot, q, q*gwf(0, ns, 1, q), xtitle = 'q (au)', title = 'Na (4p)'
oplot, q, q*gwf(0, 3.0, 1, q), linestyle = 2
oplot, q, -q*Fnl(1., 0, 3, 1, q), psym = 4
; Mg 6s
ns = 4.358
plot, q, q*gwf(0, ns, 0, q), title = 'Mg (6s)', xtitle = 'q (au)'
oplot, q, q*gwf(0, 4.0, 2, q), linestyle = 2
oplot, q, q*Fnl(1., 0., 4.0, 2, q), psym = 4
; Mg 6p
ns = 4.86
plot, q, q*gwf(0, ns, 1, q), title = 'Mg (6p)', xtitle = 'q (au)', xr=[0,1.0], /xs
oplot, q, q*gwf(0, 5.0, 1, q), linestyle = 2
oplot, q, -q*Fnl(1., 0, 5, 1, q), psym = 4
; Cs 10p
ns = 6.43
plot, q, q*gwf(0, ns, 1, q), title = 'Cs (10p)', xtitle = 'q (au)', xr=[0,1.0], /xs
oplot, q, q*gwf(0, 6.0, 1, q), linestyle = 2
oplot, q, q*Fnl(1., 0, 6, 1, q), psym = 4
print, systime(1) - t, ' SECONDS'
!p.multi = 0
!p.charsize = 1
end