;title="AR4"
;moment="ymonmean"
;fname="eir"
;step1=0
;step2=0
;countries=(/"rwanda","tanzania","kenya","uganda","burundi"/)
;tail="HF_AR4_lon27-42.0_lat-12-6.0_date19990101-19990501.nc"

;infile="/home/netapp-clima/users/tompkins/MODELS/VECTRI/VECTRI_v1.21/output/vectri_mean_rwanda_lon28-32.0_lat-4-0.0_date20001030-20101231.nc"
;infile="/home/netapp-clima/users/tompkins/MODELS/VECTRI/output_east_africa_HF/vectri_ymonstd_eafrica_lon27-42.0_lat-12-6.0_date20001030-20111001.nc"

;infile="/home/netapp-clima/users/tompkins/MODELS/VECTRI/output_east_africa_HF/vectri_mean_eafrica_lon27-42.0_lat-12-6.0_date20001030-20111001.nc"
;infile="/home/netapp-clima/users/tompkins/MODELS/VECTRI/wafrica/output/vectri_mean_wafrica_lon-18-9.0_lat4.5-19.5_date20001030-20111001.nc"
;infile="/home/netapp-clima-shared/OBS/TOPOGRAPHY/ETOPO/ETOPO2v2c_latlon.nc"
;infile="/home/netapp-clima/users/tompkins/MODELS/VECTRI/output_east_africa_HF/input/topography.nc"

infile="/home/netapp-clima/users/tompkins/MODELS/VECTRI/output_east_africa_HF/timmean_rmig0.005.nc"

outfile="./rwanda"
title=" "
countries="rwanda"
step1=-1 ; 9 
step2=-1 ; 9

; uganda
lon0=29
lon1=35
lat0=-2
lat1=5

; rwanda
lon0=28
lon1=31
lat0=-3
lat1=-1

moment="mean"
fname="pr"
