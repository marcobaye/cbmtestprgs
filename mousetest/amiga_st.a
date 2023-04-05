;ACME 0.95.7
; drivers for Amiga and Atari ST mice

; CIA #1, $dc01: %76543210
;	bit7	bit6	bit5	bit4	bit3	bit2	bi1	bit0
;Amiga:	keyb,	keyb,	keyb,	LBUT,	XB,	YB,	XA,	YA
;ST:	keyb,	keyb,	keyb,	LBUT,	YA,	YB,	XA,	XB	(swapped YA/YB names to make y axis point down)

; signal change sequences on pins 4321 (bits 3210):
; '_' means don't care, "->" shows flow of time, all sequences wrap around
;Amiga:
; up	| _0_0 -> _1_0 -> _1_1 -> _0_1
; down	| _0_0 -> _0_1 -> _1_1 -> _1_0	(y steps 0, 1, 2, 3)
; left	| 0_0_ -> 1_0_ -> 1_1_ -> 0_1_
; right	| 0_0_ -> 0_1_ -> 1_1_ -> 1_0_	(x steps 0, 1, 2, 3)
;ST:
; up	| 00__ -> 01__ -> 11__ -> 10__
; down	| 00__ -> 10__ -> 11__ -> 01__	(y steps 0, 1, 2, 3)
; left	| __00 -> __01 -> __11 -> __10
; right	| __00 -> __10 -> __11 -> __01	(x steps 0, 1, 2, 3)


	!zone amiga_st

; button poll for Amiga and Atari, called via interrupt
amiga_st_poll	; left button
		ldy #BUTTON_NOPRESS
		lda state_now
		and #$10
		bne +
			ldy #BUTTON_PRESS
+		sty sprite_amiga	; O##
		sty sprite_atarist	; O##
; CAUTION: the following only works in VICE! My real Amiga mouse connects
; potx/poty to GND, which cannot be detected by a C64. But maybe there are
; mice out there with some pullup resistor, so here's the code:
		; right button
		ldy #BUTTON_NOPRESS
		lda potx
		bpl +
			ldy #BUTTON_PRESS
+		sty sprite_amiga + 2	; ##O
		sty sprite_atarist + 2	; ##O
		; middle button
		ldy #BUTTON_NOPRESS
		lda poty
		bpl +
			ldy #BUTTON_PRESS
+		sty sprite_amiga + 1	; #O#
		sty sprite_atarist + 1	; #O#
		rts

; tables for quadrature-to-internal-format conversion
; internal format is %yycxx000
;	"c" bit is a pseudo-carry set before subtraction, so that two
;		2-bit subtractions can be done in one instruction.
;	"yy" is the step in y direction (0, 1, 2, 3 meaning 0, 1, 0, -1)
;	"xx" is the step in x direction
amiga_to_shared	!byte %........, %.#......, %....#..., %.#..#...
		!byte %##......, %#......., %##..#..., %#...#...
		!byte %...##..., %.#.##..., %...#...., %.#.#....
		!byte %##.##..., %#..##..., %##.#...., %#..#....
st_to_shared	!byte %........, %...##..., %....#..., %...#....
		!byte %##......, %##.##..., %##..#..., %##.#....
		!byte %.#......, %.#.##..., %.#..#..., %.#.#....
		!byte %#......., %#..##..., %#...#..., %#..#....

	AMIGA		= 0	; offset for tables below
	ATARI_ST	= 1
.previous	!by 0, 0	; 8 bits for Amiga, 8 for ST
.xoffs		!by ENUM_AMIGA,			ENUM_ATARIST
.yoffs		!by ENUM_AMIGA + ENUM_TOTAL,	ENUM_ATARIST + ENUM_TOTAL

amiga_st_idle ; called as often as possible
		lda state_now
		and #%....####	; keep only quadrature bits
		pha
		; Amiga
		tax
		lda amiga_to_shared, x	; convert Amiga quadrature signal to internal format
		ldy #AMIGA
		jsr .shared
		; Atari ST
		pla
		tax
		lda st_to_shared, x	; convert ST quadrature signal to internal format
		ldy #ATARI_ST
		;FALLTHROUGH
; shared part of Amiga/ST driver
; entry: A = new internal state, Y = machine index (AMIGA or ATARI_ST)
.shared		pha	; remember
		ora #%..#.....	; set pseudo carry so subtraction works as expected
		sec
		sbc .previous, y	; calculate y and x changes
; %00	no change
; %01	+1
; %10	no change (actually forbidden/invalid)
; %11	-1
		; check for y change
		asl	; C = change (inverted!), now N = "something changed!"
		bpl +
			; Y coordinate has changed
			pha
			ldx .yoffs, y
			jsr changeC
			pla
+		; check for x change
		asl	; N = pseudo-carry
		asl	; N = change (inverted!)
		asl	; C = change (inverted!), now N = "something changed!"
		bpl +
			; X coordinate has changed
			ldx .xoffs, y
			jsr changeC
+		pla
		sta .previous, y
		rts
