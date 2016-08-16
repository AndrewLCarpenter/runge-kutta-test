#!MC 1410
$!VarSet |MFBD| = '/ump/fldmd/home/mhcarpen'
$!VarSet |PATH1| = '|$PATHI|'

$!PAGE NAME = 'Untitled'
$!PAGECONTROL CREATE
$!OPENLAYOUT  "/ump/fldmd/home/mhcarpen/runge-kutta-test/tecplot scripts/err1.lay"
$!READDATASET  ' |PATH1| '
  READDATAOPTION = NEW
  RESETSTYLE = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  VARNAMELIST = '"V1" "V2"'

$!PRINTSETUP PALETTE = COLOR
$!REDRAWALL 
$!EXPORTSETUP IMAGEWIDTH = 1080
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!RemoveVar |MFBD|
