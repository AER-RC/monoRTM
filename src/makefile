SHELL = /bin/sh
CF=	f90 
LDR=	f90 

#------if we run build.sh
#COMMAND=../bin/monortm_$(vers) #if we want to automatize the version #


#------if we run make directly
vers=2.11
COMMAND=../bin/monortm_v$(vers)

SRCS=	monortm.f monortm_sub.f lblatm_monortm.f modm.f 

OBJS=	monortm.o monortm_sub.o lblatm_monortm.o modm.o 

.f.o:
	$(CF) -c $*.f

$(COMMAND): $(OBJS) 
	$(LDR) -o $(COMMAND) $(OBJS)  



