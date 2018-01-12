;*********** Load Libraries ************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$HOME/CODE/NCL/TOOLS/my_tools.ncl"

data=addfile("timmean.nc","r")
;data=addfile("1.nc","r")
;data=addfile("2.nc","r")

population=ndtooned(data->population)*1e6
eir=ndtooned(365.*data->eir)
npts=dimsizes(eir)

maxeir=500.

;print ("number "+npts+"  eir "+dimsizes(population))

;printVarSummary(eir)

bins=(/0,100,200,500,1000/)
labels=(/"0-100","100-200","200-500","500-1000","1000+"/)  

kellyhope=(/140.,105.,85.,65.,25./)
nbins=dimsizes(bins)
xbins=ispan(0,nbins-1,1)

stateir=new((/nbins,5/),"float")

stateir=0.0

pdfidx=pdfindex(population,bins)

do i=0,nbins-1
  bineir=eir(ind(pdfidx.eq.i.and.eir.lt.1000))
  stats = stat_dispersion(bineir, False)
  stateir(i,0)=stats(4)  ; LOWER OCTILE
  stateir(i,1)=stats(6)  ; LOWER QUARTILE
  stateir(i,2)=stats(0)  ; AVERAGE
  stateir(i,3)=stats(10) ; UPPER QUARTILE
  stateir(i,4)=stats(12) ; UPPER OCTILE
  delete(bineir)

;  print (stateir(i,:))
end do


;*********************************
; create plot
;********************************
 wks   = gsn_open_wks ("ps", "box" )
 res   = True
 res=sethelvetica(res)
; res@gsnScale = True               

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

; res@vpXF      = 0.10        ; In page coordinates, where to start
; res@vpYF      = 0.75        ; the plot
; res@vpHeightF = 0.43        ; Changes the aspect ratio
; res@vpWidthF  = 0.85

 res@trYMinF   =  0        ; min value on y-axis
 res@trYMaxF   =  maxeir        ; max value on y-axis

; res@tmXBMode="Explicit"
 res@tmXBLabels=labels
; res@tmXBValues=xbins
 res@tmXBLabelAngleF=45.

;**********************************************
; resources for polylines that draws the boxes
;**********************************************  
  llres                   = True			
  llres@gsLineThicknessF  = 2.5                 ; line thickness 
;**********************************************
; resources that control color and width of boxes
;**********************************************  
  opti          = True			
  opti@boxWidth = 0.7				; Width of box (x units)

 res@tiYAxisString = "Annual EIR" ; y-axis label      
 res@tiXAxisString = "Population density (km-2)" ; y-axis label      
 res@tiMainString  = "VECTRI integration for W. Africa (20km)" ; title
;***********************************************
  res@gsnMaximize = True                        ; Maximize box plot in frame.
;***********************************************

  plot=boxplot(wks,xbins,stateir,opti,res,llres)	; All 3 options used...

  mres               = True                     ; marker mods desired
  mres@gsMarkerColor = "navy blue"              ; change color
  mres@gsMarkerIndex = 1                        ; change style
  mres@gsMarkerSizeF = 80.                      ; change size


  dum=0
  do i=0,nbins-1
    str = unique_string("polymarker")
    dum@$str$ = gsn_add_polymarker(wks,plot,xbins(i),kellyhope(i),mres) 
  end do

  draw(plot)
  frame(wks)
  

; res@gsnYRefLine           = 0.              ; reference line   
; res@gsnXYBarChart         = True            ; create bar chart 
; res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
; res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

;**********************************************************************
; now, change the outline to be the color of the fill. This is useful when
; you have many data points b/c the outline will normally obscure the colors
;**********************************************************************
; res@xyLineColors = (/"red","blue"/)         ; colors for the above/blow
; plot = gsn_csm_xy (wks,xbins,aveir,res)

