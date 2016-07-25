nodes=( 0 1 2 3 4 5 6 7 8 9 10 )


for (( it=0; it<11; it++ ))
do
  SLINV="s/XVX/${nodes[it]}/"
  cat leftrightperl.template | sed $SLINV  > leftrightperl.f90
  gfortran leftrightperl.f90 -o leftrightperl
  ./leftrightperl
done
