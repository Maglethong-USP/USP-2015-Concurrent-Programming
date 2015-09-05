###############
# Directories #
###############
# Binary Output Directory
BIN			=	./bin

# Documentation Files/Directory
DOC			=	./doc

# Includes Directory
LIB			=	./lib

# Object Directiory
OBJ 		=	./obj

# Includes Directory
SRC			=	./src

#########
# Other #
#########

# Main Source File
SOURCES		=	$(SRC)/*.c

# Zip File Name
ZIP			=	project1

# Other Files to add to zip
ZIP_FILES	=	./Makefile ./README.txt ./LICENSE

# Flags
FLAGS 		=	-lm -pthread -lrt


################
# Make Modules #
################

all:
	@rm -r -v $(BIN)
	@mkdir -v $(BIN)
	gcc -I$(LIB) -o $(BIN)/prog $(SOURCES) $(FLAGS)

run:
	$(BIN)/prog

zip:
	@zip -r $(ZIP) $(LIB) $(SRC) $(DOC) $(ZIP_FILES)

clear:
	@touch $(ZIP).zip
	@rm -r $(ZIP).zip
	@touch $(BIN) $(OBJ)
	@rm -r $(BIN) $(OBJ) -v
	@touch .dummy~
	@find ./ -name *~ | xargs rm -v
