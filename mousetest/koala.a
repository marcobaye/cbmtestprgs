;ACME 0.97
; koala pad driver (doubles as paddle driver)

; CIA #1, $dc01: %76543210
;	bit7	bit6	bit5	bit4	bit3		bit2		bit1	bit0
;	keyb,	keyb,	keyb,	key,	RIGHTBUTTON,	LEFTBUTTON,	keyb,	keyb

	!zone koala

	.LEFT_BUTTON	= $04
	.RIGHT_BUTTON	= $08

koala_poll	; first x:
		ldx potx
		lda .translate_x_lo, x
		sta table_x_lo + ENUM_KOALA
		lda .translate_x_hi, x
		sta table_x_hi + ENUM_KOALA
		; now y:
		ldx poty
		lda .translate_y, x
		sta table_y_lo + ENUM_KOALA
		; left button
		lda state_now
		ldy #BUTTON_NOPRESS
		and #.LEFT_BUTTON
		bne +
			ldy #BUTTON_PRESS
+		sty sprite_koala + 0	; O##
		; right button
		lda state_now
		ldy #BUTTON_NOPRESS
		and #.RIGHT_BUTTON
		bne +
			ldy #BUTTON_PRESS
+		sty sprite_koala + 2	; ##O
		rts

; lookup tables to convert 256x256 to 297x180 (320x200 with sprites completely visible)
; these are not balanced, I was too lazy for a Bresenham...
.translate_x_lo	!for i, 0, 255 {
			!by <(i * (320-24) / 255)
}
.translate_x_hi	!for i, 0, 255 {
			!by >(i * (320-24) / 255)
}
.translate_y	!for i, 0, 255 {
			!by i * (200-21) / 255
}
