;*********** Load Libraries ************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$HOME/CODE/NCL/TOOLS/my_tools.ncl"


; will load this from plot_options
load "./plot.options.ncl"
print(infile)
vecdata=addfile(infile,"r")

jpdffile="test"

; read data
vecxa=ndtooned(vecdata->temp)
vecya=ndtooned(vecdata->rain)
vecza=ndtooned(vecdata->PRd)

goodpts=ind(vecza.ge.0.0)
vecx=vecxa(goodpts)
vecy=vecya(goodpts)

vecz=vecza(goodpts)

npts=dimsizes(vecz)
print (" number of point "+npts)


; variable 1 (e.g temperature)
xtitle="temperature"
mnx=floattointeger(floor(min(vecx)))
mxx=floattointeger(ceil(max(vecx)))
nxx=mxx-mnx+1
xbins=fspan(mnx,mxx,nxx)
nxbins=dimsizes(xbins)


; variable 2 (e.g. rain)
ytitle="rainfall"
ybins=(/0,0.1,0.2,0.5,1,2,4,8,10,50/)
nybins=dimsizes(ybins)

ztot=new((/nxbins,nybins/),"float")
ncnt=new((/nxbins,nybins/),"float")
ztot=0
ncnt=0

; indices
vecxi=pdfindex(vecx,xbins)
vecyi=pdfindex(vecy,ybins)

do i=0,npts-1
;  print (i+": "+vecxi(i)+","+vecyi(i))
  ztot(vecxi(i),vecyi(i))=ztot(vecxi(i),vecyi(i))+vecz(i)
  ncnt(vecxi(i),vecyi(i))=ncnt(vecxi(i),vecyi(i))+1.0
  avv=ztot(vecxi(i),vecyi(i))/ncnt(vecxi(i),vecyi(i))
  print ("avg "+ztot(vecxi(i),vecyi(i))+" "+ncnt(vecxi(i),vecyi(i))+" "+avv)
  if (avv.gt.1) then 
    print ("sssssssssssssssssss"+avv)
  end if 
end do

; calculate mean
do i=0,nxbins-1
do j=0,nybins-1
  if (ncnt(i,j).eq.0) then 
    ztot(i,j)=-0.1
  else
    ztot(i,j)=ztot(i,j)/ncnt(i,j)
  end if 
end do 
end do 

print (ztot)
;*********************************
; create plot
;********************************
;---Open X11 window for graphics.
  wks  = gsn_open_wks("pdf",jpdffile)
  res                      = True

; gsn_define_colormap(wks,"WhiteBlue")
;  gsn_define_colormap(wks,"GrayWhiteGray")
;  gsn_define_colormap(wks,"gsltod")

  res=mygreyscale2(wks,res)
  res=sethelvetica(res)
  res@gsnSpreadColors      = True         ;
;  res@gsnSpreadColorStart  = 127
;  res@gsnXRefLineColor="black"
;  res@gsnXRefLineColor="black"

  res@gsnMaximize          = True          ; Maximize plot in frame.
  res@cnFillOn             = True          ; Turn on contour fill
  res@cnFillMode           = "RasterFill"  ; Turn on raster fill
  res@cnLinesOn            = False         ; Turn off contour lines

;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF       = -7
;  res@cnMaxLevelValF       =  7
;  res@cnLevelSpacingF      =  1
;  res@lbLabelStride        = 2
;  res@lbBoxLinesOn         = False

  res@tmXBMode = "Explicit"
  res@tmYLMode = "Explicit"
;  res@tmXBMode = "Automatic"
;  res@tmYLMode = "Automatic"
;  res@tmXBMinorOn          = False
;  res@tmYLMinorOn          = False
  res@tmXBValues=ispan(0,nxbins-1,1)
  res@tmYLValues=ispan(0,nybins-1,1)
  res@tmXBLabels=tostring(floattointeger(xbins))
  res@tmYLLabels=(/"0","0.1","0.2","0.5","1","2","4","8","10","50"/)
print (xbins)

;  res@tmXBTickSpacingF =0.1
;  res@tmYLTickSpacingF =0.1

  res@tmXBLabelAngleF = 45
  res@tmYLLabelAngleF = 45

;  res@sfXArray = xbins
;  res@sfYArray = ybins

  res@lbOrientation="Vertical"

  res@tiXAxisString = xtitle
  res@tiYAxisString = ytitle

print (tostring(xbins))
  contour = gsn_csm_contour(wks,transpose(ztot),res)

