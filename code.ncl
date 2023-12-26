load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;--read sampling time--
ncase=17
data=asciiread("input_time.txt",-1,"string")
day=str_get_field(data(1::), 1," ")
hrsa_hhmm=stringtointeger(str_get_field(data(1::), 2," "))
hrea_hhmm=stringtointeger(str_get_field(data(1::), 3," "))
hrsc_hhmm=stringtointeger(str_get_field(data(1::), 4," "))
hrec_hhmm=stringtointeger(str_get_field(data(1::), 5," "))

;--time conversion from hhmm to float in hh
hrsa=hrsa_hhmm/100+(hrsa_hhmm%100)/60.
hrea=hrea_hhmm/100+(hrea_hhmm%100)/60.
hrsc=hrsc_hhmm/100+(hrsc_hhmm%100)/60.
hrec=hrec_hhmm/100+(hrec_hhmm%100)/60.

;--set basic info--
setfileoption ("cdf", "MissingToFillValue", False)
a=addfile("./maorwpcls."+day(0)+".cdf","r")
z=a->height
delete(a)
nz=dimsizes(z)
wa_top10=new((/ncase,nz/),float)
wa_top10@_FillValue=-9999
aerosol=new((/ncase/),float)
num_sample=new((/ncase,nz/),integer)

;--loop over all cases--
do i=0,ncase-1
  ;--read RWP data
  a=addfile("./maorwpcls."+day(i)+".cdf","r")
  timec=a->time_offset/3600
  wa=a->VerticalVelocity
  ref=a->ReflectivityUAZR
  echo=a->EchoClassification

  ;--define convection_echo_time as the time when echo classification is "convection" or "weak convection"
  nt=dimsizes(timec)
  convection_echo_time=new(nt,integer)
  do it=0,nt-1
    convection_echo_time(it)=where(any(echo(it,:).eq.6.or.echo(it,:).eq.7),1,0)
  end do
  good_time=new(nt,integer)
  do it=0,nt-1
    good_time(it)=where(any(abs(wa(it,:)).gt.40),0,1)
  end do

  ;--select updraft (1) larger than 0 m/s (2) at good time without values out of (-40,40) and (3) at convection_echo_time within convection sampling period (4) above 1.5km
  wa@_FillValue=-9999
  wa=where(wa.gt.0,wa,-9999)
ind_tc=ind(timec.ge.hrsc(i).and.timec.le.hrec(i).and.convection_echo_time.eq.1.and.good_time.eq.1)
  indz_le1_5=ind(z.le.1.5)
  wa(:,indz_le1_5)=-9999

  ;--for case 20140423, remove low level data before 1836 UTC
  if(day(i).eq."20140423") then
     indt=ind(timec.lt.18.6)
     indz_low=ind(z.le.4)
     wa(indt,indz_low)=-9999
  end if

  ;--calculate top 10th percentile averaged updraft for all heights
  do iz=0,nz-1
    if(.not.all(ismissing(wa(ind_tc,iz)))) then
      w1d_all=wa(ind_tc,iz)
      indw=ind(.not.ismissing(w1d_all))
      w1d=w1d_all(indw)
      ip=dim_pqsort(w1d, -2)
      lenw=dimsizes(w1d)
      wa_top10(i,iz)=avg(w1d(0:toint(lenw*0.1)))
      num_sample(i,iz)=dimsizes(w1d(0:toint(lenw*0.1)))
      delete([/w1d_all,w1d,indw,ip/])
    else
      wa_top10(i,iz)=-9999
      num_sample(i,iz)=0
    end if
  end do

  ;--read SMPS data
  filename=systemfunc ("ls ./maoaossmpsS1.b1."+day(i)+"*.nc")
  b=addfile(filename,"r")
  timea=b->time/3600
  dNdlogD=b->dN_dlogDp
  D=b->diameter_mobility
  D_bounds=b->diameter_mobility_bounds
  dlogD=log10(D_bounds(:,1))-log10(D_bounds(:,0))

  ;--ignore negative dNdlogD values
  dNdlogD@_FillValue=-9999
  dNdlogD=where(dNdlogD.gt.0,dNdlogD,dNdlogD@_FillValue)

  ;--calculate aerosol number concentration for D>15nm
  nta=dimsizes(timea)
  N_gt15nm=new(nta,float)
  dD=new(dimsizes(D),float)
  ind_gt15nm=ind(D.ge.15)
  do j=0,nta-1
    N_gt15nm(j)=sum(dNdlogD(j,ind_gt15nm)*dlogD(ind_gt15nm))
  end do

  ;--get aerosol concentration averaged over sampling window
  ind_ta=ind(timea.ge.hrsa(i).and.timea.le.hrea(i))
  aerosol(i)=avg(N_gt15nm(ind_ta))

  delete([/a,b,ind_tc,ind_ta,timec,timea,dNdlogD,dlogD,N_gt15nm/])
end do

;--group the top 10th averaged updrafts into 4 bins based on aerosol
aero_bin=(/0,1000,1900,3000,6000/)
nbin=dimsizes(aero_bin)-1
wa_bin=new((/nbin,nz/),float)
wa_bin@_FillValue=-9999
wa_bin!1="z"
wa_bin&z=z
num_bin=new((/nbin,nz/),integer)
se_bin=new((/nbin,nz/),float)

;--loop over bins to get averaged top 10th updrafts and corresponding standard errors for each bin
do k=0,nbin-1
  ind_bin=ind(aerosol.ge.aero_bin(k).and.aerosol.lt.aero_bin(k+1))
  wa_bin(k,:)=dim_avg_n(wa_top10(ind_bin,:),0)
  num_bin(k,:)=dim_sum_n(num_sample(ind_bin,:),0)
  ;--do smoothing for updraft plot
  indz_low=ind(z.lt.4)
  indz_up=ind(z.gt.4)
  wa_bin(k,indz_low)=runave(wa_bin(k,indz_low),9,1)
  wa_bin(k,indz_up)=runave(wa_bin(k,indz_up),17,1)

  ;--calculate standard error
  do j=0,nz-1
    if(.not.all(ismissing(wa_top10(ind_bin,j))))then
      se_bin(k,j)=stddev(wa_top10(ind_bin,j))/(sqrt(num_bin(k,j)))
    else
      se_bin(k,j)=0
    end if
  end do
  delete(ind_bin)
end do

;--do the plot
   type="x11"
  ; type="png"
   wks=gsn_open_wks(type,"wa_bin_gt15nm")
   res= True
   res@trXMinF=0
   res@trXMaxF=13
   res@trYMinF=0
   res@tiYAxisString           = "Height (km)"
   res@tiXAxisString           = "Updraft velocity (m s~S~-1 ~N~)"
   res@tiMainString      = "D > 15 nm"
   res@vpWidthF = 0.65
   res@vpHeightF = 0.45
   res@gsnStringFontHeightF    =0.025
   res@xyDashPatterns  = (/0,0,0,0/)
   res@xyLineThicknesses = (/ 12,12,12,12/)
   res@xyLineColors      = (/"black","green4","blue","red"/)
   res@gsnDraw = False
   res@gsnFrame   = False
   res@tiMainFontHeightF = 0.035
   res@tmBorderThicknessF =10
   res@tmXBMajorThicknessF = 10
   res@tmYLMajorThicknessF = 10
   res@tmXBLabelFontHeightF =0.035
   res@tmYLLabelFontHeightF =0.035
   res@gsnStringFontHeightF    =0.035
   res@tiYAxisFontHeightF =0.035
   res@tmXBMajorOutwardLengthF = 0.0
   res@tmYLMajorOutwardLengthF = 0.0
   res@tmXBMinorOn = False
   res@tmYLMinorOn = False
   plot  = gsn_csm_xy (wks,wa_bin,z,res)
   colors=(/"gray50","palegreen","LightBlue","LightPink"/)
   opacity=(/0.25,0.5,0.5,0.3/)
   dum1=new(nbin,graphic)
   dum2=new(nbin,graphic)
   gsres                   = True
   gsres@tfPolyDrawOrder   = "Predraw"
   do i=0,nbin-1
   ;--get xp,yp values for plotting shaded errors
   ind_le4km=ind(z.le.4)
   ind_gt4km=ind(z.gt.4)
   wa_bin_bk=wa_bin(i,:)
   wa_bin_bk(ind_gt4km)=-9999
   indx1=ind(.not.ismissing(wa_bin_bk))
   indx1s=min(indx1)
   indx1e=max(indx1)
   wa_bin_bk=wa_bin(i,:)
   wa_bin_bk(ind_le4km)=-9999
   indx2=ind(.not.ismissing(wa_bin_bk))
   indx2s=min(indx2)+1
   indx2e=max(indx2)-1
   nx1=indx1e-indx1s+1
   nx2=indx2e-indx2s+1
   xp1=new(nx1*2,float)
   yp1=new(nx1*2,float)
   xp2=new(nx2*2,float)
   yp2=new(nx2*2,float)
   se_bin(i,indx1s)=0
   se_bin(i,indx1e)=0
   se_bin(i,indx2s)=0
   se_bin(i,indx2e)=0
   xp1(0:nx1-1)=wa_bin(i,indx1s:indx1e)+se_bin(i,indx1s:indx1e)
   xp1(nx1::)=wa_bin(i,indx1s:indx1e)-se_bin(i,indx1s:indx1e)
   yp1(0:nx1-1)=z(indx1s:indx1e)
   yp1(nx1::)=yp1(0:nx1-1)
   xp2(0:nx2-1)=wa_bin(i,indx2s:indx2e)+se_bin(i,indx2s:indx2e)
   xp2(nx2::)=wa_bin(i,indx2s:indx2e)-se_bin(i,indx2s:indx2e)
   yp2(0:nx2-1)=z(indx2s:indx2e)
   yp2(nx2::)=yp2(0:nx2-1)
   gsres@gsFillColor       = colors(i)
   gsres@gsFillOpacityF    = opacity(i)
   dum1(i)=gsn_add_polygon (wks,plot,xp1,yp1,gsres)
   dum2(i)=gsn_add_polygon (wks,plot,xp2,yp2,gsres)
   delete([/indx1,indx2,xp1,yp1,xp2,yp2/])
   end do
   draw(plot)
   frame(wks)
end


