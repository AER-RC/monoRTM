;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;Sample Visualization software for MONORTM.
;S.A. Boukabara AER INC. 1999
;
;To run it, type: IDL>.run monortm.pro
;It will automatically look for a file
;in: out/MONORTM.OUT.
;Caution: If we run more profiles or more
;channels than this code is expecting
;('nprof' profiles and 'nchan' channels), 
;we need in this case to increase the sizes 
;of these variables before to run the idl .
;code. Sid 2001.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  PRO plot_legend,nl,xmin,xmax,ymin,ymax,linesty,col,comm,sym,sz
	  idb=0
	  nx=10
	  ny=10
	  thck=2
	  if (idb eq 0) then begin ;---linear
	   dx=(xmax-xmin)/float(nx)
	   dy=(ymax-ymin)/float(ny)
	   xleg=[xmin+(nx*0.60)*dx,xmin+(nx*0.80)*dx]
	  endif
	  if (idb eq 1) then begin ;---log
	   dx=(xmax/xmin)^(1/float(nx))
	   dy=(ymax/ymin)^(1/float(ny))
	   xleg=[xmin*dx^(nx*0.60),xmin*dx^(nx*0.80)]
	  endif
	  for i=0,nl-1 do begin
	    if (idb eq 0) then begin ;---linear
	      yleg=[ymin+(i+1)*dy,ymin+(i+1)*dy]
	    endif
	    if (idb eq 1) then begin ;---log
	      yleg=[ymin*dy^(i+1),ymin*dy^(i+1)]
	    endif
	    oplot,xleg,yleg,line=linesty(i),color=col(i),psym=-sym(i),symsize=sz,thick=thck
	    xyouts,xleg(1),yleg(1),' '+comm(i),color=col(i),charsize=sz,charthick=thck
	  endfor
	  end

;-------Main program
	loadct, 39
	device	='x'
	set_plot, 'x'
	nprof=1000	 ;we presume the maximum size (# of profiles)
	npar=18
	nchan=1000
	freq_orig	    =fltarr(nchan)
	monortm_orig	    =fltarr(npar,nchan,nprof)

;----File opening MONORTM----
	
	fic0='../out/MONORTM.OUT'
        print, 'FILE ',fic0,' OPEN FOR READ ....'
	OPENR,1,fic0
	ligne=''
	readf,1,ligne
	readf,1,ligne
	i=long(0)
	WHILE (not eof(1)) DO BEGIN
		readf,1,$
		format='(3i5,2f9.3,E19.9,f9.5,2f8.4,3f8.2,f9.3,10E12.4)',$
		IPR,NR,ich,FR,mo1,rad,mo2,W,cloud,TMPSFC,REFLC,EMISS,$
	 	ANGLE,OTOT,OTOT_WV,OTOT_O2,OTOT_N2,OTOT_O3,$
	 	OTOT_N2O,OTOT_CO,OTOT_SO2,OTOT_NO2,$
		OTOT_OH
		freq_orig(ich-1)=fr        ;frequency
		monortm_orig(0,ich-1,nr-1)=mo1      ;BRIGHTNESS TEMPERATURE (K)
		monortm_orig(1,ich-1,nr-1)=mo2      ;TOTAL TRANSMITTANCE
		monortm_orig(2,ich-1,nr-1)=cloud    ;INTEGRATED CLOUD AMOUNT (mm)
		monortm_orig(3,ich-1,nr-1)=w        ;INTEGRATED WATER VAPOR (mm)
		monortm_orig(4,ich-1,nr-1)=OTOT     ;COLUMN TOTAL OPTICAL DEPTH
		monortm_orig(5,ich-1,nr-1)=OTOT_WV  ;COLUMN WV  OPTICAL DEPTH
		monortm_orig(6,ich-1,nr-1)=OTOT_O2  ;COLUMN O2  OPTICAL DEPTH
		monortm_orig(7,ich-1,nr-1)=OTOT_N2  ;COLUMN N2  OPTICAL DEPTH
		monortm_orig(8,ich-1,nr-1)=OTOT_O3  ;COLUMN O3  OPTICAL DEPTH
		monortm_orig(9,ich-1,nr-1)=OTOT_N2O ;COLUMN N2O OPTICAL DEPTH
		monortm_orig(10,ich-1,nr-1)=OTOT_CO ;COLUMN CO  OPTICAL DEPTH
		monortm_orig(11,ich-1,nr-1)=OTOT_SO2;COLUMN SO2 OPTICAL DEPTH
		monortm_orig(12,ich-1,nr-1)=OTOT_NO2;COLUMN NO2 OPTICAL DEPTH
		monortm_orig(13,ich-1,nr-1)=OTOT_OH ;COLUMN OH  OPTICAL DEPTH
		monortm_orig(14,ich-1,nr-1)=EMISS   ;SURFACE EMISSIVITY
		monortm_orig(15,ich-1,nr-1)=TMPSFC  ;SURFACE TEMPERATURE (K)
		monortm_orig(16,ich-1,nr-1)=rad     ;RADIANCE
		monortm_orig(17,ich-1,nr-1)=ANGLE   ;ZENITH ANGLE (DEG)

		i=i+1
	ENDWHILE
	nb2=i	;number of the records read
	print, 'nb2=',nb2
	close,1
	ich0=ich-1


 sels:	
	   fl=0
	   impr=0
 	   freq=freq_orig(0:ich-1)
	   monortm=monortm_orig(*,0:ich-1,0:nr-1)

	   ;---TEST TO SEE WHICH PLOTS ARE POSSIBLE
	   IF ((ich eq 1) and (nr gt 1)) then begin
	     print, 'ONE CHANNEL PROCESSED WITH ',nr,' PROFILES'
	     tb_vs_wv      =1
	     tbspect       =0
	     radspect      =0
	     optspect_1    =0
	   ENDIF
	   ;---TEST TO SEE WHICH PLOTS ARE POSSIBLE
	   IF ((ich gt 1) ) then begin
	     print, ich,' CHANNELS PROCESSED WITH ',nr,' PROFILES'
	     tb_vs_wv      =0
	     tbspect       =1
	     radspect      =1
	     optspect_1    =1
	   ENDIF

 imp:	   
	   ;---CONTROL PARAMETERS: to select which variable to visualize.
	   ;tb_vs_wv    -->  ;to plot the tbs versus the water vapor
	   ;tbspect	-->  ;to plotthe spectrum of TB
	   ;radspect    -->  ;to plot the spectrum of radiances
	   ;optspect_1  -->  ;plots the optical depth of :h2o,o2,n2,..etc
	
	   ;---SELECT THE PROFILE TO BE PLOTTED: must be between 0 and nprof-1
	   PROFILE2PLOT=0


	   !p.font=1
	   IF (tb_vs_wv eq 1) then begin
	     erase
	     !p.multi=1
	     plot,monortm(0,0,*),monortm(3,0,*),ytitle='Water Vapor (cm)',$
	     xtitle='Brightness Temperature (K)',psym=1,symsize=0.5,$
	     charsize=1.3,charthick=4,thick=3,title='Frequency:'+string(freq(0),'(f12.3)')
	   ENDIF

	   IF (tbspect eq 1) then begin
	     erase
	     !p.multi=1
	     plot,freq,monortm(0,*,PROFILE2PLOT),xtitle='Frequency (GHz)',$
	     ytitle='Brightness Temperature (K)',$
	     charsize=1.3,charthick=4,thick=3,title='PROFILE #'+string(PROFILE2PLOT+1,'(I5)')
	   ENDIF

	   IF (radspect eq 1) then begin
	     erase
	     !p.multi=1
	     plot,freq,monortm(16,*,PROFILE2PLOT),xtitle='Frequency (GHz)',ytitle='Radiance',$
	     charsize=1.3,charthick=4,thick=3,title='PROFILE #'+string(PROFILE2PLOT+1,'(I5)')
	   ENDIF


	   IF (optspect_1 eq 1) then begin
	     erase
	     !p.multi=1
	     loadct,39
	     col_wv=140.
	     col_o2=230.
	     col_n2=90.
	     lin_wv=0
	     lin_o2=0
	     lin_n2=0
	     plot,/ylog,freq,monortm(4,*,PROFILE2PLOT),xtitle='Frequency (GHz)',$
	     ytitle='Optical Depth (Neper)',$
	     yrange=[0.0001,100],charsize=1.2,charthick=3,title='PROFILE #'+$
	     string(PROFILE2PLOT+1,'(I5)')
	     oplot,freq,monortm(5,*,PROFILE2PLOT),color=col_wv,linestyle=lin_wv
	     oplot,freq,monortm(6,*,PROFILE2PLOT),color=col_o2,linestyle=lin_o2
	     oplot,freq,monortm(7,*,PROFILE2PLOT),color=col_n2,linestyle=lin_n2
	     plot_legend,4,min(freq),max(freq),0.000001,0.006,$
	     [0,lin_wv,lin_o2,lin_n2],[0,col_wv,col_o2,col_n2],['All','H2O','O2','N2'],$
	     [0,0,0,0],0.9
	   ENDIF


 sui2:	   IF (impr eq 1) THEN BEGIN
      	   	device,/close
		device	='x'
		set_plot, 'x'
		goto, bo
	   ENDIF
	   print ,'-----------------------------------------------------'
	   boucle: print, ' PostScript File (1:YES /0:NO)'
	   read, impr
	   IF ((impr ne 1) and (impr ne 0)) THEN GOTO, boucle
	   IF (impr eq 0) THEN GOTO, bo
 ip:	   IF (impr eq 1) THEN BEGIN
		ficres='visu.ps'
		print, 'Postscript :',ficres  
		set_plot,'ps'
		device,filename=ficres,/color,ysize=23,yoffset=2,xsize=17
		goto, imp
	   ENDIF
	   bo:print ,'-----------------------------------------------------'
	   impr=0
	   device	='x'
	   set_plot, 'x' 
	   fin:close,1		;fermeture fichier
	   close,1
	   print, 'Peaceful End ..'
 	end;



