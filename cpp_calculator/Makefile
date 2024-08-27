# Makefile

SOURCE_FILE ?= $(PROJECT_ROOT)/main.cpp
CUDA_ENABLED ?= 0
PROJECT_ROOT := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))

BUILD_DIR := $(PROJECT_ROOT)build
BUILD_CLANG_DIR := $(PROJECT_ROOT)build_clang
SOURCE_REL_PATH := $(shell realpath --relative-to="$(PROJECT_ROOT)" "$(SOURCE_FILE)")
EXECUTABLE_PATH := $(BUILD_DIR)/$(SOURCE_REL_PATH:%.cpp=%)
EXECUTABLE_PATH := $(EXECUTABLE_PATH:%.cu=%)

	
ifeq ($(DEBUG),yes)
    MAKE_BUILD_TYPE = Debug
else
    MAKE_BUILD_TYPE = Release
endif

.PHONY: all
all: build

.PHONY: build
build:
	@echo $(DEBUG)
	@echo "Building project..."
	@mkdir -p $(BUILD_CLANG_DIR)
	# cmake -G Ninja -S $(PROJECT_ROOT) -B $(BUILD_CLANG_DIR) -DCMAKE_BUILD_TYPE=$(MAKE_BUILD_TYPE) -DCUDA_ENABLED=$(CUDA_ENABLED) -DUSE_CLANG=ON
	# ninja -C $(BUILD_CLANG_DIR) -t compdb > $(BUILD_CLANG_DIR)/compile_commands.json
	@mkdir -p $(BUILD_DIR)
	cmake -S $(PROJECT_ROOT) -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(MAKE_BUILD_TYPE) -DCUDA_ENABLED=$(CUDA_ENABLED) -DCMAKE_VERBOSE_MAKEFILE=0
	cmake --build $(BUILD_DIR)

.PHONY: run
run:
	@if [ -f "$(EXECUTABLE_PATH)" ]; then \
		echo "Running $(EXECUTABLE_PATH)"; \
		"$(EXECUTABLE_PATH)"; \
	else \
		echo "Executable not found: $(EXECUTABLE_PATH)"; \
	fi

.PHONY: clean
clean:
	@echo "Cleaning up..."
	@rm -rf $(BUILD_DIR)

.PHONY: debug
debug:
	@$(MAKE) build DEBUG=yes
	cp $(EXECUTABLE_PATH) $(BUILD_DIR)/debug
