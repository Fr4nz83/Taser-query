# **** GCC Toolchain path definition (better if multilib enabled...) ****
CC = gcc #C compiler;
CPP = g++ #C++ compiler
LD = g++ #Linker passthrough
OBJCOPY = objcopy # Per generare l'hex dall'elf;
OBJDUMP = objdump # Per disassemblare l'eseguibile prodotto;
NM = nm # Per generare un'analisi della dimensione delle varie sezioni degli oggetti;
ifeq (,$(OS))
RM = rm
else
RM = del
endif




# ****************** #
# Main Makefile vars #
# ****************** #
OPENCL = 0# Se "0", il progetto viene compilato escludendo il supporto OpenCL.
CUDA   = 0# Se "0", il progetto viene compilato escludendo il supporto a CUDA.
APP64  = 0# Se "1", compila l'applicativo in modalita' x64.
CUDA_INSTALL_PATH = ${CUDA_PATH}# CUDA SDK install path.
PROFILE = 0


# ************#
# Define vari #
# ************#

# OS Macro # 
ifeq (,$(OS))
DEFINES += -D_NIX_COMP_
else
DEFINES += -D_WIN_COMP_
endif

# Debug macro #
# DEFINES += -DPROJECTDEBUG




# Directories used during compilation and linking #
SRC_DIR += src# Directory sorgenti;
SRC_COMMON_DIR += srcCommon# Directory sorgenti;
BUILD_DIR += bin# Directory di linkaggio e creazione eseguibile;
LIBS_DIR += libs# Directory in cui sono presenti librerie di varia natura;
DEPS += deps# Directory dei file dipendenze;




# Main include directives #
INCLUDES = "-I$(SRC_DIR)" "-I$(LIBS_DIR)/geolib" "-I$(LIBS_DIR)/lemon"

# static/dynamic libraries directives #
#LIBDIRS = "-L$(LIBS_DIR)/geolib"
LIBS = "-lpthread" "-lboost_system" "-lboost_program_options" "-lGeographic" "-ltbb"




# *** Sources and objects to be considered by target rules ***
SRCS = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(SRC_DIR)/*.cpp)
# Tieni conto degli oggetti prodotti dal compilatore;
OBJS := $(SRCS:.c=.o)
OBJS := $(OBJS:.cpp=.o)
OBJS := $(OBJS:.cu=.o)





# *** Optimization level *** # 
# 0    = Switched off.
# 1    = Base opts;
# 2    = Default opts;
# 3    = Max opts;
# fast = Max opts with breakage of some standard ISO C++ rules;
# s    = Optimize for size;
OPT = 3
# *** x64 compile mode activation *** #
ifeq ($(APP64),1)
FLAG64 = -m64
else
FLAG64 =
endif

# *** Debug options *** #
DEBUG = 
#DEBUG += -g3

# *** Profiling options *** #
ifeq ($(PROFILE),1)
FLAGPROFILE = -pg
else
FLAGPROFILE =
endif


# **** Opzioni compilatore ****
# -c: Compila ma non linka i sorgenti; il linking viene fatto in un passo successivo sugli oggetti creati;
# -g: Genera informazioni per il debugging (vedi $(DEBUG));
# -On: Ottimizza il codice (n=0,1,2,3,s);
# -Dstringa: predefinisce "stringa" come una macro;
# -Ipath: Aggiunge un path di inclusione;
# -mcpu="fam": indica al compilatore qual'e' il target di compilazione, relativamente al codice macchina prodotto;
# -mthumb: seleziona l'instruction set adatto per il micro; l'opzione dipende dal parametro "mcpu=fam";
# -mlittle-endian: imposta la rappresentazione "little endian" in memoria per i tipi primitivi aventi dimensione >= 16 bit;
# -Wall: Attiva i warning 'normali';
# -fopenmp: attivazione supporto OpenMP (da attivare sia in fase di compilazione che di linking);
# -fsigned-char: Imposta di default la rappresentazione con segno per il tipo "char";
# -fdump-class-hierarchy: Stampa le informazioni relative al layout in memoria degli oggetti delle varie classi;
# -fdata-sections -ffunction-sections: Place each function or data item into its own section in the output file if the
#                                      target supports arbitrary sections. Permits size optimizations by the linker;
# -fomit-frame-pointer: Don't keep the frame pointer in a register for functions that don't need one (abilitato da -O);
# -ffunction-sections: ogni funzione viene collocata in una sezione specifica in .text;
# -fno-exceptions: disabilita la possibilita' di usare le eccezioni; riduce di non poco l'overhead associato alle classi;
# -fno-rtti: Disabilita la feature "Run Time Type Identifier;
# -Wl,argomento: passa l'argomento specificato al linker;
# -Wa,argomento: passa l'argomento specificato al compilatore assembler;
# -MD -MF filename.dep: per ogni sorgente, il compilatore genera le relative dipendenze e le memorizza in un file "d" apposito; utile per non fare ogni volta il build completo di un progetto;
# -MP: ai file dipendenze aggiunge anche i phony target rappresentati dagli header;
#
# *** Opzioni linker (passate via driver G++/GCC) ***
# -T"nomescript" (lato compilatore): utilizza il linker script "nomescript" indicato nell'argomento al posto di quello di default;
# --gc-sections: garbage collector del linker, rimuove i metodi/funzioni inutilizzati dall'oggetto finale;
# --print-gc-sections: stampa in stderr i metodi/funzioni "tagliati" fuori dal linker;
# 2>nomefile.txt: redireziona stderr nel file "nomefile.txt"
#
# *** Opzioni assembler ***
# -EL: genera l'oggetto tenendo conto che il processore segue l'ordine little-endian per i dati;

# Assembler, Compiler e Linker flags;
COMPILER_FLAGS = -c $(FLAGPROFILE) $(DEBUG) -O$(OPT) $(FLAG64) -std=c++11 -fopenmp -Wall -MD -MP -MF $(DEPS)/$(@F:%.o=%.d) $(DEFINES) $(INCLUDES)
CPPCOMPILER_FLAGS = -c $(FLAGPROFILE) $(DEBUG) -O$(OPT) $(FLAG64) -std=c++11 -fopenmp -Wall -MD -MP -MF $(DEPS)/$(@F:%.o=%.d) $(DEFINES) $(INCLUDES)
ASSEMBLER_FLAGS = -c $(FLAGPROFILE) $(DEBUG) -O$(OPT) $(FLAG64) -std=c++11 -fopenmp -Wall -Wa,-EL $(DEFINES)
LINKER_FLAGS = $(FLAGPROFILE) $(LIBDIRS) $(LIBS) $(FLAG64) $(DEBUG) -std=c++11 -fopenmp -Wl,--gc-sections




# **** Nome eseguibile finale ****
NOMEBIN = CASRQ
ANALISIBIN = $(NOMEBIN)Analysis
DEASMBIN = $(NOMEBIN)DeASM





# *** BEGIN CUDA check *** #
ifeq ($(CUDA),1)
OBJS := $(OBJS:.cu=.o)# genero i target dei sorgenti;
CCCUDA = nvcc

LD = nvcc
LINKER_FLAGS = $(LIBDIRS) $(LIBS) $(FLAG64)
CUDA_COMPILER_FLAGS = -c -arch=sm_20

else
DEFINES += -D_NO_CUDA_
endif
# *** END CUDA check *** #
# *** OpenCL and CUDA management *** # 




# Inclusione dei file di dipendenza (*.d) generati (eventualmente) in precedenza 
# dal compilatore.
# Queste inclusioni servono ad inserire nel makefile i target relativi ai
# componenti usati nel progetto, in maniera tale da tenere conto delle loro spe-
# cifiche dipendenze(e quindi evitare inutili ricompilazioni). 
include $(wildcard $(DEPS)/*.d)




###############
# Build project
# Major targets
###############


# target di default
# default: all


all: begin build end


begin:
	@echo Compilation start!
	@echo "Main makefile flags: CUDA $(CUDA) | OPENCL $(OPENCL) | 64BIT $(APP64)"
	@echo Translation units to compile: $(SRCS)
	@echo Static/dynamic libraries used: $(LIBS)
	@$(CC) --version
ifeq ($(CUDA),1)
ifneq (,$(OS))
		$(error CUDA can use only VS under Windows!! Error!!)
endif
endif


# Linking: creazione oggetto finale!
# in OBJS ci sono le dipendenze, ovvero gli oggetti da linkare per l'eseguibile finale;
# in LIBOBJS sono indicati gli oggetti prodotti in seguito alla compilazione della libreria FW;
# in OBJSTARTASM e' indicato l'oggetto prodotto dalla compilazione del file di startup;
build: $(LIBOBJS) $(OBJS)
# Crea l'exe;
	@echo ---------------
	@echo Linking e creazione eseguibile:
# Linka gli oggetti e crea l'eseguibile finale;
	$(LD) $^ $(LINKER_FLAGS) -o $(BUILD_DIR)/$(NOMEBIN).o
	@echo ---------------
#	@echo Conversione eseguibile nel formato eseguibile a seconda della piattaforma Win/Unix!
ifeq (,$(OS))
	$(OBJCOPY) $(BUILD_DIR)/$(NOMEBIN).o $(BUILD_DIR)/$(NOMEBIN)
else
	$(OBJCOPY) $(BUILD_DIR)/$(NOMEBIN).o $(BUILD_DIR)/$(NOMEBIN).exe
endif
# Analisi eseguibile prodotto;	
#	@echo Analisi eseguibile prodotto...
#	$(OBJDUMP) $(BUILD_DIR)/$(NOMEBIN).o -C -d > $(BUILD_DIR)/$(DEASMBIN).txt
#	$(NM) $(BUILD_DIR)/$(NOMEBIN).o -C -S --size-sort --radix=d > $(BUILD_DIR)/$(ANALISIBIN).txt

	
end:
	@echo Compilazione terminata correttamente!



######################################################################################
# Target impliciti C e C++ (usati solo se non sono state ancora generate le dipendenze
# relative alle componenti utilizzate nel progetto) 
######################################################################################
# $< indica la prima dipendenza;
# $^ indica TUTTE le dipendenze;
# $@ indica il target;
# $(@F) Indica il nome del file relativo al target (esclude quindi il path);
# $(@D) Il contrario di cui sopra; 
# Compila i sorgenti C ad uno ad uno;
#Genera le dipendenze e compila i sorgenti;
%.o: %.c
	@echo ---------------
	@echo Compilazione $<:
	$(CC) $(COMPILER_FLAGS) $< -o $@


# Compila i sorgenti CPP ad uno ad uno; 
%.o: %.cpp
	@echo ---------------
	@echo Compilazione $<:
	$(CPP) $(CPPCOMPILER_FLAGS) $< -o $@
	
	
# Compila il file di startup (assembly);
%.o: %.S
	@echo ---------------
	@echo Compilazione $<:
	$(CC) $(ASSEMBLER_FLAGS) $< -o $@
	
	
# Compila i kernel CUDA;
%.o: %.cu
	@echo ---------------
	@echo Compilazione $<:
	$(CCCUDA) $(CUDA_COMPILER_FLAGS) $< -o $@	



# Rimuovi i prodotti della compilazione e del linking;
.PHONY: clean
clean:
ifeq (,$(OS))
	@$(RM) -f $(OBJS) $(DEPS)/*.d $(BUILD_DIR)/*$(NOMEBIN)*
else
	@$(RM) $(subst /,\,$(OBJS) $(DEPS)/*.d $(BUILD_DIR)/*$(NOMEBIN)*")
endif
	@echo Cancellati oggetti, dipendenze ed eseguibili!
	
	
# Rimuovi solo i prodotti della compilazione;
cleanobjects:
ifeq (,$(OS))
	@$(RM) -f $(OBJS) $(DEPS)/*.d $(BUILD_DIR)/*$(NOMEBIN)*.o
else
	@$(RM) $(subst /,\,$(OBJS) $(DEPS)/*.d $(BUILD_DIR)/*$(NOMEBIN)*.o")
endif
	@echo Cancellati prodotti intermedi compilazione!
