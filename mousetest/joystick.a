;ACME 0.95.7
; joystick driver

; CIA #1, $dc01: %76543210
;	bit7	bit6	bit5	bit4	bit3	bit2	bi1	bit0
;	keyb,	keyb,	keyb,	BUTTON,	RIGHT,	LEFT,	DOWN,	UP

	!zone joystick

; constants
	.MAX_STEP	= $10	; ($7f max) maximum pixel step size ("speed") for joystick acceleration routine.
	.SLOW_PIXELS	= $04	; ($7f max) distance before acceleration starts, in pixels.
; vars
.pixel_counter	!byte 0	; pixel counter before accelerating


joystick_poll
		; fetch byte holding direction flags
		lda state_now
		tax	; ...and remember it
		; check 'up' direction
		ror
		bcs ++
			; subtract current step size from y value
			tay
			sec
			lda table_y_lo + ENUM_JOYSTICK
			sbc .step_size
			sta table_y_lo + ENUM_JOYSTICK
			bcs +
				dec table_y_hi + ENUM_JOYSTICK
+			tya
++		; check 'down' direction
		ror
		bcs ++
			; add current step size to y value
			tay
			;clc	; C is always clear here
			lda table_y_lo + ENUM_JOYSTICK
			adc .step_size
			sta table_y_lo + ENUM_JOYSTICK
			bcc +
				inc table_y_hi + ENUM_JOYSTICK
+			tya
++		; check 'left' direction
		ror
		bcs ++
			; subtract current step size from x value
			tay
			sec
			lda table_x_lo + ENUM_JOYSTICK
			sbc .step_size
			sta table_x_lo + ENUM_JOYSTICK
			bcs +
				dec table_x_hi + ENUM_JOYSTICK
+			tya
++		; check 'right' direction
		ror
		bcs ++
			; add current step size to x value
			tay
			;clc	; C is always clear here
			lda table_x_lo + ENUM_JOYSTICK
			adc .step_size
			sta table_x_lo + ENUM_JOYSTICK
			bcc +
				inc table_x_hi + ENUM_JOYSTICK
+			tya
++		ror
		; check button (left mouse button)
		ldy #BUTTON_NOPRESS
		bcs +
			ldy #BUTTON_PRESS
+		sty sprite_joystick	; O##
		; check "right mouse button" (for 1350 and compatibles)
		ldy #BUTTON_NOPRESS
		lda potx
		bmi +
			ldy #BUTTON_PRESS
+		sty sprite_joystick + 2	; ##O
		; check "middle mouse button" (for 1350 mode of CMD? MicroMys?)
		ldy #BUTTON_NOPRESS
		lda poty
		bmi +
			ldy #BUTTON_PRESS
+		sty sprite_joystick + 1	; #O#

; restore joystick direction bits and check whether to set speed to
; zero.
		txa
		and #$0f	; clear unneeded bits
		cmp #$0f	; any direction bit ?
		bne +
; no direction was used, so reset speed and wait counter to normal.
			lda #$01
			sta .step_size
			lda #.SLOW_PIXELS
			sta .pixel_counter
			jmp .restrict
+
; a direction bit was used, so check whether to accelerate: if speed
; is already maximum speed, don't accelerate.
.step_size = * + 1
		lda #0	; MODIFIED!
		cmp #.MAX_STEP	; if speed is max., don't accelerate
		bcs .restrict
			; speed isn't maximum yet. check whether
			; we have to wait before accelerating.
			dec .pixel_counter
			bpl .restrict
				; counter has underrun, so accelerate.
				inc .pixel_counter	; reset counter
				inc .step_size	; increase speed
.restrict
; restrict coordinate range
		ldx #ENUM_JOYSTICK
		jsr restrict
		ldx #ENUM_JOYSTICK + ENUM_TOTAL
		jmp restrict
