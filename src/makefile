CF=	f90 
LDR=	f90 
COMMAND=../bin/monortm

SRCS=	monortm.f monortm_sub.f lblatm_monortm.f modm.f 

OBJS=	monortm.o monortm_sub.o lblatm_monortm.o modm.o 

.f.o:
	$(CF) -c $*.f

$(COMMAND): $(OBJS) 
	$(LDR) -o $(COMMAND) $(OBJS)  



