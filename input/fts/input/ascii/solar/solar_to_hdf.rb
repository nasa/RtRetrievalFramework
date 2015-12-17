require 'scanf'
require 'rubygems'
require "narray"


# Simple code to write a NArray to a HDF 5 field.
def write_hdf(fname, group, data)
  # NArray is in fortran order. Rearrange to C order, and make sure this is a 
  # C double 
  data_c = data.transpose(*(0...data.dim).to_a.reverse).to_f
  File.open("t.dat", "w") do |f|
    f << data_c.to_s
  end
  system "h5import t.dat -d #{data.shape.join(",")} -p '#{group}' -t FP -s 64 -o #{fname}"
end
def write_hdf_str(fname, group, str)
  File.open("t.dat", "w") do |f|
    f << str << "\n"
  end
  system "h5import t.dat -d 1 -p '#{group}' -t STR -o #{fname}"
end

# Read in the solar data. This is made up of fixed length ASCII fields.
freqarr = [] 
strenarr = [] 
w_widarr = [] 
d_widarr = []
fname = "solar_dc_20090102.101"
File.readlines(fname).each do |ln|
# This is from the old fortran code. The data is fix length fields. We
# actually only want the first few fields. The variables aren't documented, 
# although we can guess at them.
  mw,freq,stren,w_wid,d_wid,sbhw,eprime,tdpbhw,sss =
    ln.scanf("%3c%12c%10c%10c%5c%5c%10c%8c%1c%36c")
  freqarr << freq.to_f
  strenarr << stren.to_f
  w_widarr << w_wid.to_f
  d_widarr << d_wid.to_f
end

out_fname="solar_dc_20090102.h5"
write_hdf(out_fname, 
          "Solar/Solar_Line_List/frequency", NArray.to_na(freqarr))
write_hdf(out_fname, 
          "Solar/Solar_Line_List/optical_thickness", NArray.to_na(strenarr))
write_hdf(out_fname, 
          "Solar/Solar_Line_List/folding_width", NArray.to_na(w_widarr))
write_hdf(out_fname, 
          "Solar/Solar_Line_List/doppler_width", NArray.to_na(d_widarr))
write_hdf_str(out_fname, 
              "Solar/Solar_Line_List/Source", fname)

