;*********** Load Libraries ************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$HOME/CODE/NCL/TOOLS/my_tools.ncl"


; will load this from plot_options
load "./plot.options.ncl"

print(infile)

data=addfile(infile,"r")
population=ndtooned(data->population)
eir=ndtooned(365.*data->eir)

ymax=200

; hack - will add topography to the VECTRI output.
;topodata=addfile("./input/topography_HF.nc","r")
;height=ndtooned(topodata->height)

npts=dimsizes(eir)

print ("number "+npts+"  eir "+dimsizes(population))
printVarSummary(population)
print (max(population))

bins=(/0,100,200,500,1000/)*1e-6
labels=(/"0-100","100-200","200-500","500-1000","1000+"/)  

;
; kelly hope and mckensie data:
;
kellyhope=(/140.4, 99.9 , 86.3 , 59.0 , 13.2 /)
eafrica=(/ 80.1557152821774,              37.39838808021683 ,            14.63649677367355  ,           1.953759643902479   ,  0.3598653255824266   /)     

nbins=dimsizes(bins)
xbins=ispan(0,nbins-1,1)
aveir=new(nbins,"double")
nav=new(nbins,"integer")
aveir=0.0
nav=0

pdfidx=pdfindex(population,bins)

do i=0,npts-1
  if (.not. ismissing(pdfidx(i)) .and. eir(i).lt.10000.0) then
    aveir(pdfidx(i))=aveir(pdfidx(i))+eir(i)
    nav(pdfidx(i))=nav(pdfidx(i))+1
  end if
end do

aveir(:)=aveir(:)/int2flt(nav(:))

print (aveir)
print (nav)


;*********************************
; create plot
;********************************
 wks   = gsn_open_wks ("pdf", "./eir_pop_pdf" )
 res   = True
 res = sethelvetica(res)
 res@gsnScale = True               

 ; dont advance frame
 res@gsnDraw = True
 res@gsnFrame = False

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

; res@vpXF      = 0.10        ; In page coordinates, where to start
; res@vpYF      = 0.75        ; the plot
; res@vpHeightF = 0.43        ; Changes the aspect ratio
; res@vpWidthF  = 0.85

 res@trYMinF   =  0        ; min value on y-axis
 res@trYMaxF   =  ymax        ; max value on y-axis
 res@trXMinF = -0.5    ; these are needed to prevent the plot shifting!
 res@trXMaxF = 4.5     ; 

 res@tmXBMode="Explicit"
 res@tmXBLabels=labels
 res@tmXBValues=xbins
 res@tmXBLabelAngleF=45.

 res@tiYAxisString = "Annual EIR" ; y-axis label      
 res@tiXAxisString = "Population Density (km~S~-2~N~)" ; y-axis label      
 res@tiMainString  = "" ; title

 res@gsnYRefLine           = 0.              ; reference line   

 res@gsnYAxisIrregular2Log = True

 res@gsnXYBarChart         = True            ; create bar chart 

 res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
 res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

;**********************************************************************
; now, change the outline to be the color of the fill. This is useful when
; you have many data points b/c the outline will normally obscure the colors
;**********************************************************************

 wide=0.25
 res@gsnXYBarChartBarWidth = wide

 res@gsnXYBarChartColors = (/"red"/)
 res@gsnAboveYRefLineBarPatterns=0
 res@xyLineColors = (/"black","blue"/)         ; colors for the above/blow
 plot1 = gsn_csm_xy (wks,xbins-wide,aveir,res)

 res@gsnXYBarChartColors = (/"blue"/)
 res@gsnAboveYRefLineBarPatterns=0
 res@xyLineColors = (/"black","blue"/)         ; colors for the above/blow
 plot2 = gsn_csm_xy (wks,xbins,kellyhope,res)

 res@gsnXYBarChartColors = (/"red"/)
 res@gsnAboveYRefLineBarPatterns=6
 res@xyLineColors = (/"black","blue"/)         ; colors for the above/blow
 plot3 = gsn_csm_xy (wks,xbins+wide,eafrica,res)


  lbres                    = True          ; labelbar only resources
  lbres@lbLabelFont      = "Helvetica"
  lbres@lbLabelFontHeight  = 4.0
  lbres@lbFillColors       = (/"red","red","blue"/)
  lbres@lbFillPatterns     = (/0,6,0/)
  lbres@lbPerimOn          = False
  lbres@lgPerimColor       = "white"
  lbres@vpWidthF           = 0.25     ; labelbar width
  lbres@vpHeightF          = 0.2     ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.3          ; puts space between color boxes
 barname = (/"VECTRI W.Africa","VECTRI E.Africa","Survey Data"/)
  gsn_labelbar_ndc(wks,3,barname,0.55,0.8,lbres)    ; draw right labelbar column

  frame(wks)
