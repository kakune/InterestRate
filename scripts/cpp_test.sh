cmake -S $1 -B $1/build -DCMAKE_BUILD_TYPE=Debug
cmake --build $1/build

PROJECT_ROOT="$1"
SOURCE_FILE="$2"

SOURCE_DIR=$(dirname "$SOURCE_FILE")
BUILD_DIR="$PROJECT_ROOT/build"
SOURCE_REL_PATH=$(realpath --relative-to="$PROJECT_ROOT" "$SOURCE_FILE")
EXECUTABLE_NAME=$(basename "$SOURCE_FILE" .cpp)
EXECUTABLE_PATH="$BUILD_DIR/${SOURCE_REL_PATH%.*}"

if [ -f "$EXECUTABLE_PATH" ]; then
    echo "Running $EXECUTABLE_PATH"
    "$EXECUTABLE_PATH"
else
    echo "Executable not found: $EXECUTABLE_PATH"
fi