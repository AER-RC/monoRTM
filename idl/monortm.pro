;Visualization software for MONORTM.
;S.A. Boukabara AER INC. 1999

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

	loadct, 39
;-------Main program
	device	=	'x'
	set_plot, 'x'
	taille=1 ;we presume the maximum size (# of profiles)
	npar=16
	nchan=1000
	freq_orig	    =fltarr(nchan)
	monortm_orig	    =fltarr(npar,nchan,taille)

;----File opening n1:MONORTM----
	
	fic0='../out/MONORTM.OUT'

        print, 'FILE ',fic0,' OPEN FOR READ ....'
	OPENR,1,fic0
	ligne=''
	readf,1,ligne
	readf,1,ligne
	i=long(0)
	WHILE (not eof(1)) DO BEGIN
		readf,1,$
		format='(3i5,2f9.3,f9.5,2f6.3,3f8.2,10E12.4)',$
		IPR,NR,ich,FR,mo1,mo2,W,cloud,TMPSFC,REFLC,EMISS,$
	 	OTOT,OTOT_WV,OTOT_O2,OTOT_N2,OTOT_O3,$
	 	OTOT_N2O,OTOT_CO,OTOT_SO2,OTOT_NO2,$
		OTOT_OH
		freq_orig(ich-1)=fr        ;frequency
		monortm_orig(0,ich-1,nr-1)=mo1  ;Tbs
		monortm_orig(1,ich-1,nr-1)=mo2  ;Trtot
		monortm_orig(2,ich-1,nr-1)=cloud ;cloud
		monortm_orig(3,ich-1,nr-1)=w     ;wv
		monortm_orig(4,ich-1,nr-1)=OTOT
		monortm_orig(5,ich-1,nr-1)=OTOT_WV
		monortm_orig(6,ich-1,nr-1)=OTOT_O2
		monortm_orig(7,ich-1,nr-1)=OTOT_N2
		monortm_orig(8,ich-1,nr-1)=OTOT_O3
		monortm_orig(9,ich-1,nr-1)=OTOT_N2O
		monortm_orig(10,ich-1,nr-1)=OTOT_CO
		monortm_orig(11,ich-1,nr-1)=OTOT_SO2
		monortm_orig(12,ich-1,nr-1)=OTOT_NO2
		monortm_orig(13,ich-1,nr-1)=OTOT_OH
		monortm_orig(14,ich-1,nr-1)=EMISS
		monortm_orig(15,ich-1,nr-1)=TMPSFC

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

 imp:	   
	   tbspect=1
	   optspect_1=1    ;plots for regular species:h2o,o2,n2,..etc

	   IF (tbspect eq 1) then begin
	     erase
	     !p.multi=1
	     plot,freq,monortm(0,*,0),xtitle='Frequency (GHz)',ytitle='Brightness Temperature (K)',$
	     ;charsize=1.2,charthick=3,xrange=[180,190],yrange=[160,180]
	     charsize=1.3,charthick=4,thick=3
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
	     plot,/ylog,freq,monortm(4,*,0),xtitle='Frequency (GHz)',ytitle='Optical Depth (Neper)',$
	     yrange=[0.0001,100],charsize=1.2,charthick=3
	     oplot,freq,monortm(5,*,0),color=col_wv,linestyle=lin_wv
	     oplot,freq,monortm(6,*,0),color=col_o2,linestyle=lin_o2
	     oplot,freq,monortm(7,*,0),color=col_n2,linestyle=lin_n2
	     ;plot_legend,4,min(freq),max(freq),min(monortm(4,*,0)),max(monortm(4,*,0)),$
	     ;[0,lin_wv,lin_o2,lin_n2],[0,col_wv,col_o2,col_n2],['Total','H2O','O2','N2'],[0,0,0,0],1.2
	     plot_legend,4,min(freq),max(freq),0.000001,0.006,$
	     [0,lin_wv,lin_o2,lin_n2],[0,col_wv,col_o2,col_n2],['All','H2O','O2','N2'],[0,0,0,0],0.9
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



