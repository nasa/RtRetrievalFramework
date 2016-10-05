namespace FullPhysics {
// This contains documentation for use by doxygen.
/****************************************************************//**
 \page spectrum_doxygen Spectrum related classes
 
  Note that there are a few closely related classes, with similar 
  sounding names. To keep everything straight, here is a list of
  some related classes:

  - Spectrum - This is a full spectrum, which contains a SpectralRange 
    and SpectralDomain.
  - SpectralRange - This is the 'Y' axis of the spectrum, the
    value of the spectrum at a particular point. This may have a
    jacobian and/or uncertainty value include.
  - SpectralDomain - This is the 'X' axis of the spectrum, the list of
    wavenumber/wavelength that we have values for.
  - SpectralWindow - This give the window that we use for measured
    values that we use in a retrieval (i.e., the 'low resolution'
    spectrum). This is responsible for filtering out things like bad
    pixels, point that are beyond the edges of the microwindows.
  - SpectrumSampling - This determines the high resolution points that
    we need for generating the low resolution we are fitting, i.e.,
    the spectrum before it has been convolved with the Ils. This may
    skip points if we are using something like a
    NonuniformSpectrumSampling (which is then responsible for
    interpolating to the value needed by the Ils).
  - SpectralBound - This gives the upper and lower bounds of a
    SpectralWindow. This is fairly coarse, it is used by some lower
    level classes that need to know roughly what
    wavenumber/wavelengths it will need to calculate for (e.g., LRadRt).
  - ForwardModelSpectralGrid - This contains the various grids we use
    in the ForwardModel. This has the low resolution grid (i.e., what
    passed through our SpectralWindow), the high resolution grid that
    we run the RadiativeTransfer code on (possibly with nonuniform
    spacing), and the high resolution grid interpolated to supply what
    is needed by the Ils (which might be identical the noninterpolated
    grid if we aren't doing nonuniform spacing.
*******************************************************************/
}
