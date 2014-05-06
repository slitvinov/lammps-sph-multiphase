# Here we generate the starting coordinates of the simulation
# using PACKMOL

# Running these three commands will probably take hours or even days.
# (It depends on how dense and how uniformly you need the packing to be.)
# NOTE: If PACKMOL gets stuck in an endless loop, then edit the corresponding
# "inp" file.  This should not happen.  You can also usually interrupt
# packmol after 30 minutes, and the solution at that point should be good
# enough for use.

packmol < step1_proteins.inp   # This step determines the protein's location
                               # It takes ~20 minutes (on an intel i7)
packmol < step2_innerlayer.inp # This step builds the inner monolayer
                               # It takes ~90 minutes
packmol < step3_outerlayer.inp # This step builds the outer monolayer
                               # It takes ~4 hours

# This should create a file named "system.xyz" (which moltemplate will read).
