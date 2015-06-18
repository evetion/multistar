EXTENSION = multistar
DATA = multistar--0.1.sql
PGFILEDESC = "Multistar prototype to store TINs"


MODULES = multistar
PGXS := $(shell pg_config --pgxs)
include $(PGXS)

