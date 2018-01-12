;*********** Load Libraries ************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$HOME/CODE/NCL/TOOLS/my_tools.ncl"


nb=11
xaxis=fspan(0,1,nb)
yaxis=fspan(0,1,nb)
bins=fspan(0,1,nb)
bins(nb-1)=1.e30
nbins=dimsizes(bins)-1 ; 10 bins only (1.0 included in last bin)
total=new((/nbins,nbins/),"float")


total=0

; will load this from plot_options
load "./plot.options2.ncl"

; overwrite infile
;infile="output_nlarv"+inl+"_nsurvival"+ins+"_eir_50/eir50_mean.nc"
;jpdffile="map_vectri_jpdf_nlarv"+inl+"_nsurv"+ins+".pdf"

mapfile="./input/ST_PR_mean_Eafrica.nc"

print(infile)

vecdata=addfile(infile,"r")
mapdata=addfile(mapfile,"r")

;population=ndtooned(data->population)

vecprda=ndtooned(vecdata->PRd)
mapprda=ndtooned(mapdata->pr)

goodpts=ind(vecprda.ge.0.0.and.mapprda.ge.0.0)

vecprd=vecprda(goodpts)
mapprd=mapprda(goodpts)

npts=dimsizes(vecprd)

ivec=pdfindex(vecprd,bins)
imap=pdfindex(mapprd,bins)

print (dimsizes(total)+" and "+npts)

do i=0,npts-1
;  print (i+" "+ivec(i)+" data "+vecprd(i)+mapprd(i))
  total(ivec(i),imap(i))=total(ivec(i),imap(i))+1.0
end do

total=total/sum(total)

print (sum(total))

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

  res@tmXBMode = "Manual"
  res@tmYLMode = "Manual"
  res@tmXBMinorOn          = False
  res@tmYLMinorOn          = False

  res@tmXBTickSpacingF =0.1
  res@tmYLTickSpacingF =0.1

  res@tmXBLabelAngleF = 45
  res@tmYLLabelAngleF = 45

  res@sfXArray = xaxis
  res@sfYArray = yaxis
  res@lbOrientation="Vertical"

  res@tiXAxisString = "MAP"
  res@tiYAxisString = "VECTRI"

  contour = gsn_csm_contour(wks,total,res)

