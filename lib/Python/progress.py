import sys

def write_progress(i, total, width=40, out_buffer=sys.stderr):
    # \r to return to beginning of line
    out_buffer.write('\r')
    perc_done = (i+1)/float(total)
    fmt = "[%-" + str(width) + "s] %d%%" 
    out_buffer.write(fmt % ('='*int(perc_done*width), perc_done*100))
    if (i+1) == total:
        out_buffer.write("\n")
    out_buffer.flush()
