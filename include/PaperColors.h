
#ifndef PAPERCOLORS_H
#define PAPERCOLORS_H

#include <TColor.h>

namespace PaperColors
{

    // ============================================================
    // FoCal palette
    // ============================================================

    constexpr int kFoCalBlue = 0x648EFE;
    constexpr int kFoCalViolet = 0x785EEF;
    constexpr int kFoCalMagenta = 0xDB267F;
    constexpr int kFoCalOrange = 0xFD6100;
    constexpr int kFoCalYellow = 0xFEAF02;
    constexpr int kFoCalLightYellow = 0xE9E086;
    constexpr int kFoCalGreen = 0xB7DC78;
    constexpr int kFoCalTurquoise = 0x00B6BC;
    constexpr int kFoCalLightCyan = 0x9CDAD8;

    constexpr int kFoCalPalette[] = {
        kFoCalBlue,
        kFoCalViolet,
        kFoCalMagenta,
        kFoCalOrange,
        kFoCalYellow,
        kFoCalLightYellow,
        kFoCalGreen,
        kFoCalTurquoise,
        kFoCalLightCyan};

    // ============================================================
    // Run comparison
    // ============================================================

    constexpr int kRun3 = 0xFF0000; // Red
    constexpr int kRun2 = 0x654BEB; // Blue

    // ============================================================
    // Particle species
    // ============================================================

    constexpr int kKstar = kFoCalOrange;  // Orange
    constexpr int kPhi = kFoCalTurquoise; // Teal

    // ============================================================
    // Monte Carlo models
    // Palette B: Pythia blue / EPOS purple-teal / Thermal orange
    // ============================================================

    constexpr int kMonash = 0x16A1C0;    
    constexpr int kShoving = 0x1A63BE;   
    constexpr int kRopes = 0x648EFE;     
    constexpr int kEPOSUrQMD = 0x8EB000; 
    constexpr int kEPOS = 0x35B000;      
    constexpr int kThermal = 0xE87522;   

    // ============================================================
    // xT scaling: collision energies
    // Sunset palette
    // Nile Blue -> Pine -> Sea Green -> Tangerine -> Orange
    // ============================================================

    constexpr int kXT2760 = kFoCalLightCyan;
    constexpr int kXT5020 = kFoCalTurquoise;
    constexpr int kXT7000 = kFoCalGreen;
    constexpr int kXT8000 = kFoCalLightYellow;
    constexpr int kXT13000 = kFoCalYellow;
    constexpr int kXT13600 = kFoCalOrange;

    // ============================================================
    // xT exponent: energy ratios
    // Saturated sunset / turquoise variation
    // ============================================================

    constexpr int kXTRatio1 = kFoCalOrange;
    constexpr int kXTRatio2 = kFoCalTurquoise;
    constexpr int kXTRatio3 = kFoCalYellow;
    constexpr int kXTRatio4 = kFoCalBlue;
    constexpr int kXTRatio5 = kFoCalYellow;
    constexpr int kXTRatio6 = kFoCalLightYellow;

    // ============================================================
    // ROOT color conversion
    // ============================================================

    inline int GetROOTColor(int hex)
    {
        return TColor::GetColor(
            (hex >> 16) & 0xFF,
            (hex >> 8) & 0xFF,
            hex & 0xFF);
    }

    // ============================================================
    // Convenience arrays
    // ============================================================

    constexpr int kXTEnergies[] = {
        kXT2760,
        kXT5020,
        kXT7000,
        kXT8000,
        kXT13000,
        kXT13600};

    constexpr int kXTRatios[] = {
        kXTRatio1,
        kXTRatio2,
        kXTRatio3,
        kXTRatio4,
        kXTRatio5,
        kXTRatio6};

    constexpr int kModels[] = {
        kMonash,
        kShoving,
        kRopes,
        kEPOSUrQMD,
        kEPOS,
        kThermal};

}

#endif
