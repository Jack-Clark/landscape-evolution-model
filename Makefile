NVCC          := nvcc

MODULES       := FA_kernels FD_kernels MEM_kernels MOD_kernels .
SRC_DIR       := $(MODULES)
BUILD_DIR     := $(addprefix build/,$(MODULES))

SRC           := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cu))
OBJ           := $(patsubst %.cu,build/%.o,$(SRC))
HEADERS       := headers $(CUDA_HOME)/include $(CUDA_HOME)/samples/common/inc
INCLUDES      := $(addprefix -I,$(HEADERS))


build/%.o: %.cu
	$(NVCC) $(INCLUDES) -c $< -o $@

.PHONY: all checkdirs clean

all: checkdirs build/lem

build/lem: $(OBJ)
	$(NVCC) $^ -o $@ -lgdal


checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf build
