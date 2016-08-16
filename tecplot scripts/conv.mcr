#!MC 1410
$!VarSet |MFBD| = '/ump/fldmd/home/mhcarpen/runge-kutta-test/TMP'
$!VarSet |PATH1| = '|$PATHI|'
$!VarSet |PATH2| = '|$PATHE|'
$!PAGE NAME = 'Untitled'
$!PAGECONTROL CREATE

$!IF |$VARS| == 2
$!OPENLAYOUT "/ump/fldmd/home/mhcarpen/runge-kutta-test/tecplot scripts/conv2.lay"
$!ELSEIF |$VARS| == 3
$!OPENLAYOUT "/ump/fldmd/home/mhcarpen/runge-kutta-test/tecplot scripts/conv3.lay"
$!ENDIF

$!READDATASET  ' |PATH1| |PATH2| '
  READDATAOPTION = NEW
  RESETSTYLE = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  VARNAMELIST = '"V1" "V2"'
$!PICK ADDATPOSITION
  X = 2.37511436414
  Y = 3.20013723696
  CONSIDERSTYLE = YES
$!REDRAWALL 
$!EXPORTSETUP IMAGEWIDTH = 1080
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!RemoveVar |MFBD|
#!RemoveVar |PATH1|
#!RemoveVar |PATH2|
