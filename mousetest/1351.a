;ACME 0.95.7
; 1351 driver (and compatibles, which might have a middle button)

; CIA #1, $dc01: %76543210
;	bit7	bit6	bit5	bit4		bit3	bit2	bi1		bit0
;	keyb,	keyb,	keyb,	LEFTBUTTON,	keyb,	keyb,	MIDDLEBUTTON,	RIGHTBUTTON

	!zone 1351

	.LEFT_BUTTON	= $10
	;.MIDDLE_BUTTON	= $02	; only for CMD mice and MicroMys
	;.RIGHT_BUTTON	= $01
	.X_STUFF	= 0
	.Y_STUFF	= 1

; variables
; (not really words: first byte is x value, second byte is y value)
.minimum	!word $ffff	; lowest pot values yet
.limit		!word 0	; highest pot values yet, plus 1
.previous	!word 0	; old values
.width		!word 0	; interval width
.half_width	!word 0	; (buffered for speed increase)

cbm1351_poll
		; mouse x
		ldx #.X_STUFF
		jsr .pot_delta
		; now YYAA is signed x movement. add to current x value
		clc
		adc table_x_lo + ENUM_1351
		sta table_x_lo + ENUM_1351
		tya
		adc table_x_hi + ENUM_1351
		sta table_x_hi + ENUM_1351
		ldx #ENUM_1351
		jsr restrict
		; mouse y
		ldx #.Y_STUFF
		jsr .pot_delta
		; now YYAA is signed y movement. subtract from current y
		; value because axis points the other way.
		clc
		sbc table_y_lo + ENUM_1351
		eor #$ff
		sta table_y_lo + ENUM_1351
		tya
		sbc table_y_hi + ENUM_1351
		eor #$ff
		sta table_y_hi + ENUM_1351
		ldx #ENUM_1351 + ENUM_TOTAL	; "y stuff"
		jsr restrict
	; mouse buttons
		; left mouse button
		lda state_now
		ldy #BUTTON_NOPRESS
		and #.LEFT_BUTTON
		bne +
			ldy #BUTTON_PRESS
+		sty sprite_1351	+ 0	; O##
		; right mouse button
		lda state_now
		lsr	; right button is in lsb (low-active)
		ldy #BUTTON_NOPRESS
		bcs +
			ldy #BUTTON_PRESS
+		sty sprite_1351 + 2	; ##O
		; middle mouse button (CMD and MicroMys only!)
		lsr	; right button is in "down" bit (low-active)
		ldy #BUTTON_NOPRESS
		bcs +
			ldy #BUTTON_PRESS
+		sty sprite_1351 + 1	; #O#
		rts

; on entry: X is direction handle (0 = x, 1 = y)
; on exit: YYAA is signed distance
.pot_delta ; compute signed distance of mouse movement
		; first, get new value and clear "recalculate signal width" flag
		lda potx, x
		ldy #0
		; check whether new value is lower than lowest known
		cmp .minimum, x
		bcs +
			; store new "lowest" und set "recalculate signal width" flag
			sta .minimum, x
			dey;ldy#$ff
+		; check whether new value is higher than highest known
		cmp .limit, x
		bcc ++
			; set "recalculate signal width" flag and store new "highest"
			ldy #$ff
			pha	; remember
			adc #0	; add one (C set)
			sta .limit, x
			; value $ff (0 after adding) means that there is no mouse connected,
			; so reset min/max in that case.
			bne +
				; no mouse, so reset "lowest"
				; ("highest" will have been reset already)
				; and return zero.
				tay		; set Y to zero.
				pla		; fix stack
				lda #$ff	; reset "lowest"
				sta .minimum, x
				tya		; return YYAA = 0
				rts
+			pla	; restore
++		; if flag is set, recalculate signal width
		iny	; check flag
		bne ++
			tay		; remember
			lda .limit, x	; get highest + 1
			sec		; subtract lowest
			sbc .minimum, x
			bcc +
				sta .width, x		; store signal
				lsr			; width and half signal
				sta .half_width, x	; width
+			tya	; restore
++		; calculate distance
		tay	; remember
		sec
		sbc .previous, x
		pha
		tya
		sta .previous, x
		pla
		beq .zero	; if not moved, exit
		bcc .negative	; negative difference
;.positive
		; check whether movement caused a value wrap-around
		cmp .half_width, x
		bcc .decrease
		beq .decrease
		; it did, so calculate "real" distance and jump to exit
		;sec	; C is always set here
		sbc .width, x	; fix distance
		; we now know that the (fixed) distance is really negative, so
		; we finally wipe out that annoying bit 0 noise by incrementing
		; the value.
.increase	;clc	; C is always clear here
		adc #1
		beq .zero	; if increasing gives zero, jump to zero handler
		ldy #$ff	; set up high byte for negative values
		rts
.negative
		; check whether movement caused a value wrap-around
		eor #$ff	; complement
		cmp .half_width, x
		eor #$ff	; restore value
		bcc .increase
		; it did, so calculate "real" distance and exit
		clc
		adc .width, x	; fix distance
		; we now know that the (fixed) distance is really positive, so
		; we finally wipe out that annoying bit 0 noise by decrementing
		; the value
.decrease	sec
		sbc #1
		; no difference or positive difference; both need zero as the high byte.
.zero		ldy #0
		rts
