Run chmod command to give run permissioin for two scripts first, then run two scripts.
compile_adjoint.sh will report an error of missing some libs, mx will fix it and genertate "geos". Move "geos" to your geoschem run directory.

chmod +x compile_adjoint.sh mx
./compile_adjoint.sh
./mx
