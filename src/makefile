CF=	f90 
LDR=	f90 
#CF=	f77 
#LDR=	f77 
COMMAND=../bin/monortm

SRCS=	monortm.f monortm_sub.f lblatm.f modm.f 

OBJS=	monortm.o monortm_sub.o lblatm.o modm.o 

.f.o:
	$(CF) -c -static $*.f

$(COMMAND): $(OBJS) 
	$(LDR) -o $(COMMAND) $(OBJS)  



