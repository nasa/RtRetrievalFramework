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

[ "kahn", "merra"].each do |prop_dir|
  hdf_tmp = "aerosol_uncompress_" + prop_dir + ".h5"
  Dir.glob(prop_dir + "/*.mie") do |f|
    name = File.basename(f, ".mie")
    print prop_dir + "/" + "#{name}.mie" << "\n"
  # Read all the lines 
    t = File.readlines(prop_dir + "/" + "#{name}.mie")
  # skipping the first one line
    t.shift
  # Go through the remaining lines to get wavenumber, qext and qsca
    wn = []
    qext = []
    qsca = []
    t.reverse!.each do |ln|
      wl, qextv, qscav = ln.split
      if(wl)
        wn << 1e4 / wl.to_f
        qext << qextv.to_f
        qsca << qscav.to_f
      end
    end

    # Write to HDF file
    write_hdf(hdf_tmp, "#{name}/Properties/wave_number", NArray.to_na(wn))
    write_hdf(hdf_tmp, "#{name}/Properties/extinction_coefficient", 
              NArray.to_na(qext))
    write_hdf(hdf_tmp, "#{name}/Properties/scattering_coefficient", 
              NArray.to_na(qsca))
    
    # Read the phase function moments
    pflist = []
    File.open(prop_dir + "/" + "#{name}.mom") do |f|
      f.readline                  # Skip first line
      until(f.eof?)
        ln = f.readline
        wl, num = ln.split
        pf = []
        (num.to_i + 1).times { pf << (f.readline.split.map {|d| d.to_f}) }
        pf_siewert = NArray.to_na(pf).transpose(1,0)
        raise "Need to have 6 columns" unless(pf_siewert.shape[1] ==6)
        # Convert to de Rooij convention
        pf = NArray.float(*pf_siewert.shape)
        pf[true,0] = pf_siewert[true,2-1]
        pf[true,1] = pf_siewert[true,1-1]
        pf[true,2] = pf_siewert[true,6-1]
        pf[true,3] = pf_siewert[true,4-1]
        pf[true,4] = -pf_siewert[true,3-1]
        pf[true,5] = pf_siewert[true,5-1]
        pf[0,0] = 1.0
        pflist << pf
      end
    end
    pflist.reverse!               # Change from wavelength to wavenumber order
    raise "Phase functions should be same size as mie" unless(pflist.size == wn.size)
    nrow = (pflist.map {|p| p.shape[0]}).max
    pfmom = NArray.float(pflist.size, nrow, 6)
    pfmom[true,true,true] = 0
    pflist.each_with_index do |pf, i|
      pfmom[i,0...(pf.shape[0]), true] = pf[true,true]
    end
    # Write to HDF file
    write_hdf(hdf_tmp, "#{name}/Properties/phase_function_moment", pfmom)
    write_hdf_str(hdf_tmp, 
                  "#{name}/Properties/Phase Function Moment Convention", "de Rooij")
    write_hdf_str(hdf_tmp, 
                  "#{name}/Properties/Source", name)
  end

  # Finally, compress the file
  system "h5repack -i " + hdf_tmp + " -o aerosol_" + prop_dir + ".h5 -f GZIP=9"
end

