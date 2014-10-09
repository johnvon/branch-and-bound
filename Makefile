SYSTEM  = x86-64_sles10_4.1

LIBFORMAT = static_pic

####diretorios com parser
TSPPARSERDIR  = ../tsp-parser
METAHDIR      = ../metaheuristics

#### define o compilador
CPPC = g++
#############################

#### opcoes de compilacao e includes
CCOPT      = -std=c++0x -g -fPIC -fexceptions -DIL_STD -Wall -O3
TSPINCLUDE = $(TSPPARSERDIR)/include
TSPLIB     = $(TSPPARSERDIR)/obj

METAHINCLUDE = $(METAHDIR)/include
METAHLIB     = $(METAHDIR)/obj

CCFLAGS    = $(CCOPT) -I$(TSPINCLUDE) -I$(METAHINCLUDE) -L$(METAHLIB) -L$(TSPLIB)
#############################

#### flags do linker
CCLNFLAGS = -std=c++0x -g -I$(TSPINCLUDE) -I$(METAHINCLUDE) -L$(METAHLIB) -L$(TSPLIB) -ltsp_parser -lmetah -O3
#############################

#### diretorios com os source files e com os objs files
SRCDIR = src
OBJDIR = obj
INCDIR = include
#############################

#### lista de todos os srcs e todos os objs
SRCS = $(wildcard $(SRCDIR)/*.cpp)
HDRS = $(wildcard $(INCDIR)/*.h)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
#############################

#### regra principal, gera o executavel
bnb: $(OBJS) 
	@echo "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(HDRS) $(OBJS) $(CCLNFLAGS) -o $@ 
############################

#inclui os arquivos de dependencias
-include $(OBJS:.o=.d)

#regra para cada arquivo objeto: compila e gera o arquivo de dependencias do arquivo objeto
#cada arquivo objeto depende do .c e dos headers (informacao dos header esta no arquivo de dependencias gerado pelo compiler)
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC) -std=c++0x -MM -I$(TSPINCLUDE) -I$(METAHINCLUDE) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

#delete objetos e arquivos de dependencia
clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d

rebuild: clean bnb
