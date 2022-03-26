#csh

setenv CMDDIR ~/CMDfit
setenv CMDOBJ  ${CMDDIR}/Obj
setenv CMDMOD  ${CMDDIR}/Mod
setenv CMDEXE  ${CMDDIR}/bin
setenv CMDDATA ${CMDDIR}/Data

# Use -i4 and -r4 on 64 bit machines.
# The if statement protects a startup of the ark from being overwritten for FC and FCC.
if ($?FC == 0) setenv FC  "/usr/local/bin/gfortran -fno-second-underscore -O2 -I$CMDMOD"
if ($?FCC == 0) setenv FCC "/usr/local/bin/gfortran -ggdb -fcheck=all -fno-second-underscore -fbounds-check -fbacktrace -I$CMDMOD"
setenv G95_ENDIAN LITTLE
setenv ARK_FITS_BSWAP TRUE

alias binup     $CMDEXE/binup  
alias subset    $CMDEXE/subset  
alias uncer     $CMDEXE/uncer
alias for_gaia  $CMDEXE/for_gaia  
alias grid      $CMDEXE/grid  
alias monte     $CMDEXE/monte  
alias tau2      $CMDEXE/tau2
alias old_tau2  $CMDEXE/old_tau2
alias sim       $CMDEXE/sim
alias iso       $CMDEXE/iso
alias ascii2cluster $CMDEXE/ascii2cluster
