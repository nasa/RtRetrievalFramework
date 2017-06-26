; based on screenread and saveimage procedures in CH 8 of "Practical
; IDL Programming" by Gumley
pro write_png_file, filename
    tvrd_true = !d.flags and 128
    if (tvrd_true eq 0) then $
      message, 'TVRD not supported on this device: ' + !d.name

    device, get_visual_depth=depth, get_decomposed=entry_decomposed

    if (depth le 8) then begin
        image = tvrd(true=0)
        tvlct, r, g, b, /get
        write_png, filename, image, r, g, b
    endif else begin
        device, decomposed=1
        image = tvrd(true=1)
        device, decomposed=entry_decomposed
        write_png, filename, image
    endelse
end
