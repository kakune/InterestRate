cmake -S $1 -B $1/build -DCMAKE_BUILD_TYPE=Release
cmake --build $1/build
$1/build/src/main