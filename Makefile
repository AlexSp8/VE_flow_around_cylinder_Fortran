# Compiler and flags
FC = mpiifort

FFLAGS = -O0 -g -traceback -check all -check bounds -check uninit -ftrapuv \
				 -gen-interfaces -warn interfaces -debug all -implicitnone \
				 -assume realloc_lhs -fstack-protector -assume protect_parens

# FFLAGS = -O0 -g
FFLAGS = -O3
# FFLAGS = -O3 -traceback -static -fp-model strict

# PETSc and SLEPc includes and libs
PETSC_INC = -I/home/aspyridakis/My_Libraries/petsc_intel/include
PETSC_LIB = -L/home/aspyridakis/My_Libraries/petsc_intel/lib -lpetsc -lm -qopenmp

SLEPC_INC = -I/home/aspyridakis/My_Libraries/slepc_intel/include
SLEPC_LIB = -L/home/aspyridakis/My_Libraries/slepc_intel/lib -lslepc -lm -qopenmp

LIBS = $(PETSC_LIB) $(SLEPC_LIB)

EXE = exe_alex
BUILD_DIR = build

# List all sources here, or read from a file if you want
SRC = \
	src/Fem3D_RCM.F90 \
	src/Parameters_mod.F90 \
	src/Global_Variables_mod.F90 \
	src/Tools_mod.F90 \
	src/Storage_mod.F90 \
	src/Constitutive_Models_mod.F90 \
	src/Constitutive_Models_SQRT_mod.F90 \
	src/Physical_mod.F90 \
	src/FEM_mod.F90 \
	src/Boundary_mod.F90 \
	src/GaussIntegration/Gauss_Line_mod.F90 \
	src/GaussIntegration/Gauss_Triangle_mod.F90 \
	src/GaussIntegration/Gauss_Quadrangle_mod.F90 \
	src/GaussIntegration/Gauss_Tetrahedron_mod.F90 \
	src/GaussIntegration/Gauss_Hexahedron_mod.F90 \
	src/GaussIntegration/Gauss_mod.F90 \
	src/Petsc_mod.F90 \
	src/Unknowns_Arrays_mod.F90 \
	src/Element_Calculations_mod.F90 \
	src/Residuals_mod.F90 \
	src/Equations_mod.F90 \
	src/Extra_Equations_mod.F90 \
	src/Dirichlet_mod.F90 \
	src/Newton_Raphson_mod.F90 \
	src/Projection_mod.F90 \
	src/Post_Process_mod.F90 \
	src/Stability/Constitutive_Models_Stability_mod.F90 \
	src/Stability/Constitutive_Models_SQRT_Stability_mod.F90 \
	src/Stability/Element_Calculations_Stability_mod.F90 \
	src/Stability/Residuals_Stability_mod.F90 \
	src/Stability/Equations_Stability_mod.F90 \
	src/Stability/Extra_Equations_Stability_mod.F90 \
	src/Stability/Post_Process_Stability_mod.F90 \
	src/Stability/Stability_mod.F90 \
	src/Stability/Continuation_Stability_mod.F90 \
	src/Continuation_Transient_mod.F90 \
	src/Continuation_Arclength_mod.F90 \
	src/Continuation_mod.F90 \
	src/Main.F90

# Generate object file names (mirror folder structure in build/)
OBJ = $(patsubst src/%.F90,$(BUILD_DIR)/%.o,$(SRC))
OBJ := $(patsubst src/Stability/%.F90,$(BUILD_DIR)/Stability/%.o,$(OBJ))
OBJ := $(patsubst src/GaussIntegration/%.F90,$(BUILD_DIR)/GaussIntegration/%.o,$(OBJ))

INC = $(PETSC_INC) $(SLEPC_INC)

# Python dependency generator script location
DEPS_SCRIPT = fortran_deps.py
SRC_LIST = src_files.txt
DEPS_FILE = deps.mk

# Default target
all: $(EXE)

# Link executable
$(EXE): $(OBJ)
	@echo $(FFLAGS)
	@echo "Linking executable $@"
	@$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Compile source to object
$(BUILD_DIR)/%.o: src/%.F90
	@echo "Compiling $< → $@"
	@mkdir -p $(dir $@)
	@$(FC) $(FFLAGS) $(INC) -c $< -o $@ -module $(BUILD_DIR)

$(BUILD_DIR)/Stability/%.o: src/Stability/%.F90
	@echo "Compiling $< → $@"
	@mkdir -p $(dir $@)
	@$(FC) $(FFLAGS) $(INC) -c $< -o $@ -module $(BUILD_DIR)

$(BUILD_DIR)/GaussIntegration/%.o: src/GaussIntegration/%.F90
	@echo "Compiling $< → $@"
	@mkdir -p $(dir $@)
	@$(FC) $(FFLAGS) $(INC) -c $< -o $@ -module $(BUILD_DIR)

# Clean build files and executable
clean:
	rm -rf $(BUILD_DIR) $(EXE) $(DEPS_FILE) $(SRC_LIST)

.PHONY: all clean deps

# Write src_files.txt (list of sources) for dependency script
$(SRC_LIST):
	@echo "Generating source list file..."
	@printf "%s\n" $(SRC) > $(SRC_LIST)

# Generate dependency file using Python script
$(DEPS_FILE): $(SRC_LIST) $(DEPS_SCRIPT)
	@echo "Generating Fortran module dependencies..."
	@python3 $(DEPS_SCRIPT) $(SRC_LIST) $(BUILD_DIR) > $(DEPS_FILE) || (rm -f $(DEPS_FILE); exit 1)

# Include generated dependencies (ignore if missing)
-include $(DEPS_FILE)

# Target to regenerate dependencies manually
deps: $(DEPS_FILE)
	@echo "Dependencies regenerated."
