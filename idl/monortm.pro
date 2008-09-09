;PRO monortm, monortm_orig, freq_orig, npar, nchan, nprof, nmol, outfil=outfil

; Modified: August 2008 VHP Update for new Monortm v4.0 output format
; Modified: 10-FEB-06  VHP Changed formatting of reading to match changes made in MONORTM 
;                           
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;Sample Visualization software for MONORTM.
;S.A. Boukabara AER INC. 1999
;
;To run it, type: IDL>.run monortm.pro
;It will automatically look for a file
;in: ../out/MONORTM.OUT.
;Caution: If we run more profiles or more
;channels than this code is expecting
;('nprof' profiles and 'nchan' channels), 
;we need in this case to increase the sizes 
;of these variables before running the idl .
;code. Sid 2001.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



;-------Main program
;	loadct, 39
	device	='x'
	set_plot, 'x'
	nprof=100	 ;we presume the maximum size (# of profiles)
	npar=9+39+1     ; 9 parameters + 39 ods by molecule + xsec od
	nchan=14000
	freq_orig	    =fltarr(nchan)
	monortm_orig	    =fltarr(npar,nchan,nprof)

;----File opening MONORTM----
	
	IF NOT( KEYWORD_SET(OUTFIL)) THEN fic0='../out/MONORTM.OUT' $
          ELSE fic0 = outfil
        
        print, 'FILE ',fic0,' OPEN FOR READ ....'
	OPENR,1,fic0
	ligne=''
	readf,1,ligne
        readf, 1,ligne
	readf,1,format='(5x,I8)',nwn
        readf, 1, ligne

        test =  strsplit(ligne, ' ')
        np =  n_elements(test)
        nmol =  np-13
	i=long(0)
	WHILE (not eof(1)) DO BEGIN

		readf,1,ligne
                dum =  fltarr(np)
                reads, ligne, dum
                ipr = dum(0)
		freq_orig(i)=dum(1)        ;frequency
		monortm_orig(0,i,ipr-1)=dum(2)      ;BRIGHTNESS TEMPERATURE (K)
		monortm_orig(1,i,ipr-1)=dum(3)     ;RADIANCE
		monortm_orig(2,i,ipr-1)=dum(4)      ;TOTAL TRANSMITTANCE
		monortm_orig(3,i,ipr-1)=dum(5)        ;INTEGRATED WATER VAPOR (cm)
		monortm_orig(4,i,ipr-1)=dum(6)    ;INTEGRATED CLOUD AMOUNT (mm)
		monortm_orig(5,i,ipr-1)=dum(7)  ;SURFACE TEMPERATURE (K)
		monortm_orig(6,i,ipr-1)=dum(8)   ;SURFACE EMISSIVITY

		monortm_orig(7,i,ipr-1)=dum(10)   ;ZENITH ANGLE (DEG)
		monortm_orig(8,i,ipr-1)=dum(11)     ;COLUMN TOTAL OPTICAL DEPTH

                FOR idum= 12,12+nmol-1 DO $
                  monortm_orig(idum-3, i, ipr-1)=dum(idum) ; ods by molecule  
		monortm_orig(12+nmol-3,i,ipr-1)= dum(12+nmol)  ;column xsec od

		i=i+1
	ENDWHILE
	nb2=i	;number of the records read
;	print, 'nb2=',nb2
	close,1

        gf =  where(freq_orig NE 0., ngf)

        set_plot, 'x'
        !p.multi =  0
        !p.region =  0
        plot, freq_orig(gf), monortm_orig(0,gf,0), $
          xtitle= 'frequency [GHz]', $
          ytitle= 'Brightness temperature [K]', $
          psym=1

 	end;



