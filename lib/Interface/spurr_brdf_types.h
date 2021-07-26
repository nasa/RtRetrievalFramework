#ifndef SPURR_BRDF_TYPES_H
#define SPURR_BRDF_TYPES_H

namespace FullPhysics {

    /// These match numbers used internal to the Fortran code for Spurr based RTs
    enum SpurrBrdfType {
        LAMBERTIAN  = 1,
        ROSSTHIN    = 2,
        ROSSTHICK   = 3,
        LISPARSE    = 4,
        LIDENSE     = 5,
        HAPKE       = 6,
        ROUJEAN     = 7,
        RAHMAN      = 8,
        COXMUNK     = 9,
        BREONSOIL   = 10,
        BREONVEG    = 11,
        BPDFNDVI    = 12,
        NEWCMGLINT  = 13,
        RTKHOTSPOT  = 14,
        MODFRESNEL  = 15,
        SNOWBRDF    = 16
    };
 
}

#endif
