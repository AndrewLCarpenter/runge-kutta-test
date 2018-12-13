 #! /bin/sh
PATHTOP='/ump/fldmd/home/mhcarpen/runge-kutta-test/Code/'

PATHIMEX=$PATHTOP'IMEX/'
PATHIMPLICIT=$PATHTOP'IMPLICIT/'

#CASES=('ESDIRK436L2SA_ARK' 'QESDIRK436L2SA' 'ESDIRK436L2SA' 'ESDIRK43I6L2SA' 'ESDIRK536L2SA' 'ESDIRK436L2SA_C3' 'ARK436[2]SA_1' 'ARK436[2]SA_2' 'ARK436[2]SA_3' 'ARK436[2]SA_4' 'ARK436[1]SA_1' 'ARK436[1]SA_2' 'ARK436[1]SA_3' 'Lirk4' 'ARK437[2]SA_1-C2EC2I' 'ARK537[2]SA_1' 'ARK548[2]SA_1' 'ARK548[2]SA_2')
#PROBLEMS=('vanderPol' 'Pureschi' 'Kaps' 'Kreiss' 'Lorenz' 'Rossler_3' 'Brusselat' 'Burgers' 'Boscar_31' 'Broadwell')
#VAR_LIST=(     2           2       2        2        3         3          2          1         2            3      )
CASES=('ESDIRK324L[2]SA' 'ESDIRK325L[2]SA' 'SDIRK4' 'ESDIRK436L[2]SA_1' 'ESDIRK436L[2]SA_2' 'ESDIRK437L[2]SA' 'ESDIRK547L[2]SA_1' 'ESDIRK547L[2]SA_2' 'ESDIRK548L[2]SA')
PROBLEMS=('vanderPol' 'Pureschi' 'Kaps' 'Kreiss' 'Lorenz' 'Boscar_31' 'Broadwell')
VAR_LIST=(     2           2       2        2        3        3           3      )
for cas in "${CASES[@]}"; do
  for prob in ${!PROBLEMS[@]}; do
    PATHFILE=${PROBLEMS[$prob]}'/'$cas'/'${PROBLEMS[$prob]}'_'$cas'_conv.dat'
    export VARS=${VAR_LIST[$prob]}
#   export PATHI=$PATHIMEX$PATHFILE
    export PATHI=$PATHIMPLICIT$PATHFILE
    export PATHE=$PATHIMPLICIT$PATHFILE    
    tec360 -b -p conv.mcr
    PLOTNAME='Plots/'${PROBLEMS[$prob]}'/'${PROBLEMS[$prob]}'_'$cas'_conv.png'
    mv export.out $PLOTNAME
    
    for ((ii=1; ii<=$VARS; ii++)); do
      PATHFILE=${PROBLEMS[$prob]}'/'$cas'/'${PROBLEMS[$prob]}'_'$cas'_'$ii'.dat'
      export PATHI=$PATHIMEX$PATHFILE
      tec360 -b -p err.mcr
      PLOTNAME='Plots/'${PROBLEMS[$prob]}'/'${PROBLEMS[$prob]}'_'$cas'_'$ii'.png'   
      mv export.out $PLOTNAME
    done

  done
done
