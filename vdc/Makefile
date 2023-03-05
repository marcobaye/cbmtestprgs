PROGS		= vdc-interactive.prg
PETCAT		= petcat
FLAGS		= -w70 -l 1c01
RM		= rm

all: $(PROGS)

vdc-interactive.prg: interactive.bas.txt
	$(PETCAT) $(FLAGS) -o vdc-interactive.prg interactive.bas.txt

clean:
	-$(RM) -f *.o *.tmp $(PROGS) *~ _*.tmp* core
