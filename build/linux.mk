# Shell Tools for Linux

# File seperator
SEPARATOR := /

# Seperator for PATH variable
PATH_SEPARATOR := :

# Create a file or overwrite if it exists
TOUCH := touch

# Copy files
COPY := cp

# Rename files
define RENAME =
$(shell mv $(1) $(2) 2>/dev/null)
endef

# Delete files
DEL := rm -f

# Look recursively for directories
define LOOK_DIR =
$(shell find $(1) -type d -ls)
endef

# Create directory only if it does not exist
define MKDIR =
$(shell mkdir -p $(1))
endef

# Copy direcotries
CPDIR := cp -r

# Delete directory in a recursive way
RMDIR := rm -r -f
