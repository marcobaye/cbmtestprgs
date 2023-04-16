ASSEMBLER6502	= acme
AS_FLAGS	= -v2 -Wtype-mismatch -r report.a
RM		= rm

PROGS		= mousetest64.prg mousetest128.prg
SRCS		= 1351.a amiga_st.a cx22.a joystick.a koala.a mousetest.a sprites.a

all: $(PROGS)

mousetest64.prg: mt64.a $(SRCS)
	$(ASSEMBLER6502) $(AS_FLAGS) --outfile mousetest64.prg --format cbm mt64.a

mousetest128.prg: mt128.a $(SRCS)
	$(ASSEMBLER6502) $(AS_FLAGS) --outfile mousetest128.prg --format cbm mt128.a

clean:
	-$(RM) -f *.o *.tmp $(PROGS) *~ core
