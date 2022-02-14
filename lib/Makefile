F2PY = f2py
F2PY_FLAGS = --fcompiler=intelem --compiler=intelem

OBJ_read_file = read_file.f90

read_file: $(OBJ_read_file)
	$(F2PY) -c $(OBJ_read_file) -m read_file $(F2PY_FLAGS)