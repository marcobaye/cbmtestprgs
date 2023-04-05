ASSEMBLER6502	= acme
AS_FLAGS	= -v2 -Wtype-mismatch -r report.a
RM		= rm

PROGS		= mousetest.prg
SRCS		= 1351.a amiga_st.a cx22.a joystick.a mousetest.a sprites.a

all: $(PROGS)

mousetest.prg: $(SRCS)
	$(ASSEMBLER6502) $(AS_FLAGS) --outfile mousetest.prg --format cbm mousetest.a

clean:
	-$(RM) -f *.o *.tmp $(PROGS) *~ core
