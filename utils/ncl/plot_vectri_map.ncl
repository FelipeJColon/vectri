;*********** Load Libraries ************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$HOME/CODE/NCL/TOOLS/my_tools.ncl"

begin 

; default options
latlonlab=15.0
lonname="longitude"
latname="latitude"
diff=0
ftitle=" "

; will load this from plot_options
load "./plot.options.ncl"

if (fname.eq."pr") then
;  varname="cases_ratio"
  varname="PRd"
  fsf=1.
  ftitle="Parasite Ratio"
  units=""
  flevels=(/0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
end if

if (fname.eq."eir") then
  varname="eir"
  fsf=365.
  ftitle="EIR"
  units="bites per year"
  flevels=(/0.1,1,5,10,20,30,50,70,100,150,200,300,500/)
end if

if (fname.eq."lts") then
  varname="lts"
  fsf=1.
  ftitle="lts"
  units="days"
  flevels=(/10,20,30,60,90,120,150,180,210,240,270,300,330,365/)
end if

if (fname.eq."cs" .or. fname.eq."cs1"  .or. fname.eq."cs1bin") then
  varname="PRd" ; hasnt been cleaned yet
  fsf=1.
  ftitle="cs"
  units=""
  flevels=(/0.1,0.3,0.5,0.7,0.9,1.0/)
end if

if (fname.eq."cs2" .or. fname .eq. "cs2bin") then
  varname="PRd"
  fsf=1.
  ftitle="cs class"
  units=""
  flevels=(/0.1,0.3,0.5,0.7,0.9,1.0/)*2
end if

if (fname.eq."height") then
  varname="height"
  fsf=1.
  ftitle="Height"  
  units="m"
;  flevels=(/250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000/)
  flevels=fspan(500,2200,18)
end if
if (fname .eq. "mappr") then 
  varname="pr"
  fsf=1.
  ftitle="Parasite Ratio"
  units=""
  flevels=(/0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
end if 
if (fname .eq. "rain") then 
  varname="rain"
  fsf=1.
  ftitle="Rainfall"
  units="mm/day"
;  flevels=tofloat((/0.1,0.5,1,1.5,2,3,5,7,10,15,20,30/))
  flevels=tofloat((/0.1,0.3,0.5,1,1.5,2,3,5,7,10/))
end if 

if (fname .eq. "temp" .or. fname .eq. "tas" .or. fname .eq. "t2m") then 
  varname="t2m"
  fsf=1.
  ftitle="Temperature"
  units="oC"
  flevels=fspan(15,30,16)
end if 

ftitle0=ftitle

; use half the range, subtract the minimum and reverse if difference plot
; NOTE - flevels must have an even size!
reverse_colmap=0
if (diff.eq.1) then 
  flevh=flevels(0:dimsizes(flevels)/2-1)
  flevels=array_append_record(-flevh(::-1),flevh,0)
  colmap="nrl_sirkes"
  if (fname .eq. "rain") then 
    reverse_colmap=1
  end if
  ; special treatment for temperature  
  if (fname .eq. "temp" .or. fname .eq. "tas" .or. fname .eq. "t2m") then 
    delete(flevels)
    flevels=fspan(0,3,31)
    colmap="WhiteYellowOrangeRed"
  end if 
else
  colmap="WhiteBlueGreenYellowRed"
end if 
if (moment.eq."ymonstd".or.moment.eq."std") then
  flevels=flevels/5
end if

;print (infile)
data=addfile(infile,"r")
lon=data->$lonname$
lat=data->$latname$
;if (fname.eq."height" .or. fname.eq."mappr") then 
;  time=2010
;else
;  lon=data->longitude
;  lat=data->latitude
  if (fname.eq."height" .or. fname.eq."mappr") then 
    time=2010
  else 
    time=data->time
  end if
;end if

if (lat(0).eq.max(lat)) then 
  latreverse=-1
else 
  latreverse=1
end if 
lat=lat(::latreverse)

if (lon0.eq.-99) then 
  delete(lon0)
  lon0=min(lon)
end if
if (lon1.eq.-99) then 
  delete(lon1)
  lon1=max(lon)
end if
if (lat0.eq.-99) then 
  delete(lat0)
  lat0=min(lat)
end if
if (lat1.eq.-99) then
  delete(lat1)
  lat1=max(lat)
end if

wks = gsn_open_wks("pdf",outfile)          ; open a ncgm file
;gsn_define_colormap(wks,"gui_default")
gsn_define_colormap(wks,colmap)

do step=step1,step2
  if (step.eq.0) then
    fld=dim_avg_n_Wrap(data->$varname$,0)
    datestr=" "
  else
    if (step.ge.1) then 
      fld3=data->$varname$
      fld=fld3(step-1,:,:)
      delete(fld3)
      itime=time(step-1)
      sdate=new(3,"integer")
      sdate=datesep(itime)
      delete(itime)
      datestr=month_str(sdate(1))
    else 
      fld=data->$varname$
    end if
  end if 
  fld=fld(::latreverse,:)

  if (title.eq."default") then 
  if (fname.eq."mappr") then
    title=title+" - "+ftitle
  else
    if (step.eq.0) then 
      title=title+" "+moment+" "+ftitle
    else
      title=ftitle0+" "+moment+" "+datestr
    end if
  end if
  end if

;print (fld)

  fld=fld*fsf
  fld@long_name=title
  fld@units=units
  clevels=flevels

  maxfld=max(fld)
  minfld=min(fld)

  res                   = True
  res=sethelvetica(res)
  res@cnRasterModeOn         = True
  res@cnFillOn               = True               ; color fill
  res@cnLinesOn              = True               ; no contour lines
  res@cnLineLabelsOn         = True               ; no contour labels
  res@cnLineLabelInterval    = 1
  res@cnInfoLabelOn          = False              ; no contour info label
  res@mpGridAndLimbDrawOrder    = "PreDraw"
  res@cnLabelDrawOrder          = "PreDraw"

  res@gsnRightStringOrthogonalPosF=0.0
  res@gsnRightStringParallelPosF = 1.2

  res@cnFillDrawOrder = "Predraw" ; important for masking
  res@cnLineDrawOrder = "Predraw" ; important for masking


  res@gsnSpreadColors        = True               ; use total colormap
;  res@gsnSpreadColorStart    = 1
  res@gsnSpreadColorEnd      = -1
  res@gsnAddCyclic           = False              ; regional data

  if (reverse_colmap.eq.1) then 
    gsn_reverse_colormap(wks)
  end if 
  res@pmTickMarkDisplayMode  = "Always"           ; turn on tickmarks


; !!!!! any plot of data that is on a native grid, must use the "corners"
; method of zooming in on map.


  res@cnLevelSelectionMode="ExplicitLEvels"
  res@cnMinLevelValF=1.
  res@cnLevels = clevels

  ; vertical bar
  res@lbOrientation          = "vertical" 

  res@gsnMaximize = True                        ; Maximize box plot in frame.

;************************************************
; mp resources
; Turn on masking, and set resource indicating
; areas to mask.
;************************************************
  res@mpAreaMaskingOn       = 1
  res@mpMaskAreaSpecifiers  = countries
;************************************************
; Turn off the fill boundaries, since we are going
; to explicitly set which boundaries we want to
; fill.
;************************************************
  res@mpFillBoundarySets    = "NoBoundaries"
;  fill_specs = (/"water","land"/)
; res@mpFillAreaSpecifiers  = fill_specs

  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  =  lat0
  res@mpLeftCornerLonF  =  lon0
  res@mpRightCornerLatF =  lat1
  res@mpRightCornerLonF =  lon1

  res@mpPerimOn           = True
  res@mpOutlineBoundarySets = "National"
  res@mpGridAndLimbOn        = True               ; lat/lon lines on
  res@mpGridSpacingF         = latlonlab               ; lat/lon lines every 10 deg
  res@mpGridLineDashPattern  = 2                  ; lat/lon lines as dashed

  res@mpDataResolution="Finest"
  res@mpDataBaseVersion = "MediumRes"            ; choose higher resolution
  res@mpDataSetName = "Earth..4"                 ; choose most recent boundaries
  res@tmXBMode="Manual"
  res@tmYBMode="Manual"
  res@tmXBTickSpacingF=latlonlab
  res@tmYBTickSpacingF=latlonlab
;
; title
;
;  res@tiMainString       =title   ; add a title
;  res@tiMainFontHeightF  = .018                               ; font height
;print (lat0)

contour = gsn_csm_contour_map(wks,fld,res)  ; create the plot
;print ("here")

end do

end
