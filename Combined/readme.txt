Combinedv_1 =   Total Absorptions * Energy Stored, for various radii over height range
                Uses F = 25nN, coat = 10nm
                Max height when force gives bending of 1/2 pitch
                Depends on Bending_Limits, Light_Absorptions_Test, Function Test

combinedv_2 =   Same as above, but 'normalised' absorption within TiO2 only.
                Depends on Function_Test, and BL2, LAT2 (any 2 labelled objects)

combinedv_2_rom =   Same as v2, with normalised absorption
                    Uses Rule-Of-Mixtures to better calculate bending/energies
                    Uses transverse mixing formula for matrix/fibre approximation

combinedv_3_rom =   Same as above, but comparing pitches to height.

combinedv_4_rom =   Same as above, but comparing pitches and radius.


function_test = Has 1 function - calculate_e(rad, height)
                Calculates energy for given radius, height under force 25nN
                (Has in-built Youngs modulus calculator.)
                Note - does not (currently) include effects of coating. This is more of a guide.

bending_limits= Func1 returns rs, heights, energies
                Calculates max height allowed for each radius (+ energy at max height)

LAT         =   func2 Calculates absorption for coat = 10, pitch=200, height/rad specified.
                Return dict with k = radius of pillar + coat, v = arr([heights tested (100-hmax), corresponding Abs))