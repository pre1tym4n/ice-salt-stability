;
; Iteration 1
; 

spjyr=31557600d0

zz=[0,25d-6,9d-5,0.0001,0.0004,0.002,0.003,0.005,0.007,0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
zz=[0,0.002,0.003,0.005,0.007,0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
pp=dblarr(zz.length)
tt=dblarr(zz.length)
pp1=dblarr(zz.length)
tt1=dblarr(zz.length)
pp2=dblarr(zz.length)

wsalt=0.10d

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp}) & $
  pp1[i+1]=q.pp & $
  tt[i+1]=q.tt[-1] & $
endforeach

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp1}) & $
  pp2[i+1]=q.pp & $
  tt1[i+1]=q.tt[-1] & $
endforeach

diff_010=(tt1-tt)



wsalt=0.20d

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp}) & $
  pp1[i+1]=q.pp & $
  tt[i+1]=q.tt[-1] & $
endforeach

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp1}) & $
  pp2[i+1]=q.pp & $
  tt1[i+1]=q.tt[-1] & $
endforeach

diff_020=(tt1-tt)



wsalt=0.3d

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp}) & $
  pp1[i+1]=q.pp & $
  tt[i+1]=q.tt[-1] & $
endforeach

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp1}) & $
  pp2[i+1]=q.pp & $
  tt1[i+1]=q.tt[-1] & $
endforeach

diff_030=(tt1-tt)


wsalt=0.5d

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp}) & $
  pp1[i+1]=q.pp & $
  tt[i+1]=q.tt[-1] & $
endforeach

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp1}) & $
  pp2[i+1]=q.pp & $
  tt1[i+1]=q.tt[-1] & $
endforeach

diff_050=(tt1-tt)


wsalt=0.75d

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp}) & $
  pp1[i+1]=q.pp & $
  tt[i+1]=q.tt[-1] & $
endforeach

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp1}) & $
  pp2[i+1]=q.pp & $
  tt1[i+1]=q.tt[-1] & $
endforeach

diff_075=(tt1-tt)


wsalt=0.9d

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp}) & $
  pp1[i+1]=q.pp & $
  tt[i+1]=q.tt[-1] & $
endforeach

foreach e, zz[1:-1], i do begin & $
  q=salty_ice_loss1(e,wsalt=wsalt,pp_profile={z:zz,pp:pp1}) & $
  pp2[i+1]=q.pp & $
  tt1[i+1]=q.tt[-1] & $
endforeach

diff_090=(tt1-tt)

id=[0,lindgen(16-6+1)+6]
w=window(dimensions=[800,800])
p=objarr(10)
p[0]=plot(zz[id],diff_010[id]/spjyr/1d6,axis_style=1,clip=0,/current,name='10 wt%')
p[0].symbol='o' & p[0].sym_filled=1
p[0].xthick=2 & p[0].xtickdir=1 & p[0].xticklen=0.012 & p[0].xtickfont_name='Arial' & p[0].xtickfont_size=18 & p[0].xtitle='Depth (m)'
p[0].ythick=2 & p[0].ytickdir=1 & p[0].yticklen=0.012 & p[0].ytickfont_name='Arial' & p[0].ytickfont_size=18 & p[0].ytitle='Delay (Myr)'
p[1]=plot(zz[id],diff_020[id]/spjyr/1d6,clip=0,symbol='s',sym_filled=1,name='20 wt%',/over)
p[2]=plot(zz[id],diff_050[id]/spjyr/1d6,clip=0,symbol='tu',sym_filled=1,name='50 wt%',/over)
p[3]=plot(zz[id],diff_075[id]/spjyr/1d6,clip=0,symbol='tl',sym_filled=1,name='75 wt%',/over)
p[4]=plot(zz[id],diff_090[id]/spjyr/1d6,clip=0,symbol='D',sym_filled=1,name='90 wt%',/over)
p[0].sym_size=0.5
p[1].sym_size=0.5
p[0].yrange=[0,60]
ll=legend(target=p,position=[0.15,65],/data)
ll.font_size=16
ll.transparency=100
p[0].position=[0.11,0.09,0.95,0.95]

w.background_color='black'
p[0].color='white'
p[0].xcolor='white'
p[0].ycolor='white'
p[1].color='white'
p[2].color='white'
p[3].color='white'
p[4].color='white'
ll.color='white'
ll.text_color='white'
p[0].xtickfont_size=24
p[0].ytickfont_size=24
p[0].position=[0.14,0.11,0.95,0.95]


end